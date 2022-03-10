import logging
import numpy as np
import pandas as pd
from itertools import chain, combinations
import matplotlib

matplotlib.use('Agg')

CNV_COLUMNS = ['caller', 'chromosome', 'start', 'end', 'status', 'cn', 'cn1',
               'cn2']
CNV_DTYPES = {'caller': str, 'chromosome': str, 'start': np.float64,
              'end': np.float64, 'status': str, 'cn': np.float64,
              'cn1': np.float64, 'cn2': np.float64}


def init_logging(appname, level=logging.INFO, output_filename=None):
    """
    Initialize console logger.

    :param appname: Name of the app for which logging is initialized
    :param level: Log-level, default INFO
    :param output_filename: If provided, logging is sent to console AND to file
    :return: Initialized Logger object
    """

    # create logger
    if output_filename is not None:
        logging.basicConfig(filename=output_filename, level=level)
    logger = logging.getLogger(appname)
    logger.setLevel(level)

    # # create file handler which logs even debug messages
    # fh = logging.FileHandler('appname.log')
    # fh.setLevel(logging.DEBUG)

    # create console handler with INFO log level
    ch = logging.StreamHandler()
    ch.setLevel(level)

    # create formatter and add it to the handler(s)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    # add handler(s) to the logger
    # logger.addHandler(fh)
    logger.addHandler(ch)

    return logger


class Region:
    """
    Container for region data. Can be derived from a single truth or query row,
    or describe overlap of any number of truth and query rows.

    sources: keeps a list of source rows for each caller represented in this
             region
    """

    def __init__(self, chromosome, start, end, sourcess,
                 allow_duplicates=False):
        """
        Initialize a new Region object.

        :param chromosome: Chromosome for region
        :param start: Start of region, inclusive
        :param end: End of region, exclusive
        :param sourcess: List of dicts of references to source rows
        :param allow_duplicates: Duplicate callers allowed in the region
        """
        self.chromosome = chromosome
        self.start = start
        self.end = end
        # This merges input list of dicts sourcess into self.sources dict
        self.sources = {}
        # Run through list of dicts
        for sources in sourcess:
            for caller in sources:
                if not allow_duplicates and caller in self.sources:
                    msg = (
                        'Cannot add duplicate {} to region {}({:010d}, '
                        '{:010d}).'.format(caller, chromosome,
                                           int(start), int(end)
                                           )
                    )
                    raise UserWarning(msg)

                self.sources[caller] = sources[caller]

    def __repr__(self):
        """
        Region is represented by the following format (start and end are left-
        padded with 0's to maintain 10 digit integers):
        chromosome(start, end)[callers present in region].

        Notes:
        - Caller listing is alphabetically sorted.
        - CAUTION: This function is used as a key for sorting regions, be
          careful when modifying!

        Example:
        chrX(0000012345, 0000067890)['BED', 'CNVkit', 'PureCN', 'truth']

        :return: String representation of this Region
        """
        return '{}({:010d}, {:010d}){}'.format(self.chromosome,
                                               int(self.start),
                                               int(self.end),
                                               sorted(self.sources.keys()))

    def intersect(self, region, allow_duplicates=False):
        """
        Intersect two regions to determine head, overlap and tail regions.

        CAUTION: Assumes 'self' comes before 'region',
                 sorting should be performed outside!

        :param region: Region to intersect with self region
        :param allow_duplicates: Duplicate callers allowed in the region
        :return: Outs list containing head and overlap regions (may be an empty
                 list), tail region (may be None)
        """
        outs = []
        tail = None

        # If different chromosomes, then self region is head, no overlap and
        # other region is tail (assumes 'self' comes before 'region', sorting
        # should be performed outside!)
        if self.chromosome != region.chromosome:
            outs.append(self)
            tail = region
            return outs, tail

        # Head
        if self.start < region.start:
            outs.append(
                Region(self.chromosome, self.start,
                       min(self.end, region.start),
                       [self.sources]))
        elif self.start > region.start:
            outs.append(Region(region.chromosome, region.start,
                               min(region.end, self.start), [region.sources]))

        # Overlap
        if self.end > region.start and self.start < region.end:
            outs.append(Region(self.chromosome, max(self.start, region.start),
                               min(self.end, region.end),
                               [self.sources, region.sources],
                               allow_duplicates))

        # Tail
        if self.end < region.end:
            tail = Region(region.chromosome, max(self.end, region.start),
                          region.end, [region.sources])
        elif self.end > region.end:
            tail = Region(self.chromosome, max(self.start, region.end),
                          self.end, [self.sources])

        return outs, tail

    def to_series(self, callers):
        """
        Get a Pandas Series representation of this region to be used to form
        the analysis DataFrame.

        :param callers: List of callers to take into account
        :return: Pandas Series representation of this region
        """
        data = {'chromosome': self.chromosome,
                'start': self.start,
                'end': self.end,
                'length': self.end - self.start}

        for caller in sorted(callers):
            # Set caller status column
            status = np.NaN
            if caller in self.sources:
                status = self.sources[caller]['status']
            data['{}_status'.format(caller)] = status

        return pd.Series(data)

    def to_series_with_cn(self, callers):
        """
        Get a Pandas Series representation of this region to be used to form
        the analysis DataFrame.

        :param callers: List of callers to take into account
        :return: Pandas Series representation of this region
        """
        data = {'chromosome': self.chromosome,
                'start': self.start,
                'end': self.end,
                'length': self.end - self.start}

        for caller in sorted(callers):
            # Set caller cn column
            cn_status = np.NaN
            if caller in self.sources:
                cn_status = self.sources[caller]['cn']
            data['{}_cn_status'.format(caller)] = cn_status

        return pd.Series(data)

    def to_series_with_cr(self, callers_and_ploidy_dict):
        """
        Get a Pandas Series representation of this region to be used to form
        the analysis DataFrame.

        :param callers_and_ploidy_dict: Dict of callers and ploidy values
        :return: Pandas Series representation of this region
        """
        callers = callers_and_ploidy_dict.keys()
        data = {'chromosome': self.chromosome,
                'start': self.start,
                'end': self.end,
                'length': self.end - self.start}

        for caller in sorted(callers):
            ploidy = callers_and_ploidy_dict[caller]
            # Set caller cn column
            cr_status = np.NaN
            if caller in self.sources:
                cr_status = round((self.sources[caller]['cn'] / ploidy), 2)
            data['{}_cr_status'.format(caller)] = cr_status

        return pd.Series(data)

    def to_series_for_consensus(self):
        """
        Get a Pandas Series representation of this region to be used to form
        the consensus DataFrame.

        Notes:
        - 'freec': caller name for Control-FREEC input file
        - 'gatk_called': caller name for GATK4 CNV called.seg input file
        - 'gatk_model': caller name for GATK4 CNV modelFinal.seg input file

        :return: Pandas Series representation of this region
        """
        data = {'chromosome': self.chromosome, 'start': int(self.start),
                'end': int(self.end)}

        # Init Control-FREEC columns
        freec_status = np.NaN
        freec_genotype = np.NaN
        freec_cn = np.NaN
        freec_cn1 = np.NaN
        freec_cn2 = np.NaN
        freec_prediction = 'neutral'

        if 'freec' in self.sources:
            freec_status = self.sources['freec']['status']
            freec_genotype = self.sources['freec']['genotype']
            freec_cn = self.sources['freec']['cn']
            freec_cn1 = self.sources['freec']['cn1']
            freec_cn2 = self.sources['freec']['cn2']
            freec_prediction = freec_status

        # Init GATK called columns
        gatk_status = np.NaN
        gatk_cn = np.NaN
        gatk_call = np.NaN
        gatk_ml2cr = np.NaN
        gatk_prediction = 'neutral'

        if 'gatk_called' in self.sources:
            gatk_status = self.sources['gatk_called']['status']
            gatk_cn = self.sources['gatk_called']['cn']
            gatk_call = self.sources['gatk_called']['CALL']
            gatk_ml2cr = self.sources['gatk_called']['MEAN_LOG2_COPY_RATIO']
            gatk_prediction = gatk_status

        # Init GATK model columns
        gatk_mafp50 = np.NaN

        if 'gatk_model' in self.sources:
            gatk_mafp50 = self.sources['gatk_model'][
                'MINOR_ALLELE_FRACTION_POSTERIOR_50']

        # Make a (consensus) prediction
        if freec_prediction == gatk_prediction:
            status = freec_prediction
            confidence = 'high'
            support = 'freec,gatk'
        elif freec_prediction == 'neutral':
            status = gatk_prediction
            confidence = 'low'
            support = 'gatk'
        elif gatk_prediction == 'neutral':
            status = freec_prediction
            confidence = 'low'
            support = 'freec'
        else:
            status = 'complex'
            confidence = 'no'
            support = np.NaN

        # Prepare output dict
        data['status'] = status
        data['confidence'] = confidence
        data['support'] = support

        data['freec_status'] = freec_status
        data['freec_genotype'] = freec_genotype
        data['freec_cn'] = freec_cn if np.isnan(freec_cn) else str(
            int(freec_cn))
        data['freec_cn1'] = freec_cn1 if np.isnan(freec_cn1) else str(
            int(freec_cn1))
        data['freec_cn2'] = freec_cn2 if np.isnan(freec_cn2) else str(
            int(freec_cn2))

        data['gatk_status'] = gatk_status
        data['gatk_cn'] = gatk_cn
        data['gatk_call'] = gatk_call
        data['gatk_mean_log2_copy_ratio'] = gatk_ml2cr
        data['gatk_minor_allele_fraction_posterior_50'] = gatk_mafp50

        return pd.Series(data)

    def to_series_for_gatk(self):
        """
        Get a Pandas Series representation of this region to be used to form
        the GATK DataFrame (combination of called.seg and modelFinal.seg files)

        Notes:
        - 'gatk_called': caller name for GATK4 CNV called.seg input file
        - 'gatk_model': caller name for GATK4 CNV modelFinal.seg input file

        :return: Pandas Series representation of this region
        """
        data = {'chromosome': self.chromosome,
                'start': self.start,
                'end': self.end,
                'length': self.end - self.start}

        # Init GATK called columns
        status = np.NaN
        cn = np.NaN
        ml2cr = np.NaN
        call = np.NaN

        if 'gatk_called' in self.sources:
            status = self.sources['gatk_called']['status']
            cn = self.sources['gatk_called']['cn']
            ml2cr = self.sources['gatk_called']['MEAN_LOG2_COPY_RATIO']
            call = self.sources['gatk_called']['CALL']

        # Init GATK model columns
        mafp50 = np.NaN
        l2crp50 = np.NaN

        if 'gatk_model' in self.sources:
            mafp50 = self.sources['gatk_model'][
                'MINOR_ALLELE_FRACTION_POSTERIOR_50']
            l2crp50 = self.sources['gatk_model'][
                'LOG2_COPY_RATIO_POSTERIOR_50']

        # Infer cn1 and cn2 (LOH)
        cn1 = np.NaN
        cn2 = np.NaN

        if not np.isnan(cn) and not np.isnan(mafp50):
            cn1 = cn * mafp50
            cn2 = cn - cn1

        data['status'] = status
        data['cn'] = cn if np.isnan(cn) else round(cn, 1)
        data['cn1'] = cn1 if np.isnan(cn1) else round(cn1, 1)
        data['cn2'] = cn2 if np.isnan(cn2) else round(cn, 1) - round(cn1, 1)
        data['call'] = call
        data['mean_log2_copy_ratio'] = ml2cr
        data['log2_copy_ratio_posterior_50'] = l2crp50
        data['minor_allele_fraction_posterior50'] = mafp50

        return pd.Series(data)

    def to_series_for_sv(self):
        """
        Get a Pandas Series representation of this region to be used to form
        the SV initiated DataFrame.

        CAUTION: region.sources contains SVs 'DEL' and 'DUP', not caller names!

        :return: Pandas Series representation of this region
        """
        data = {'chromosome': self.chromosome,
                'start': self.start,
                'end': self.end}

        if len(self.sources.keys()) == 1:
            if 'DEL' in self.sources:
                # Assume minimal loss
                data['status'] = 'loss'
                data['cn'] = 1
                data['cn1'] = 0
                data['cn2'] = 1
            if 'DUP' in self.sources:
                # Assume minimal gain
                data['status'] = 'gain'
                data['cn'] = 3
                data['cn1'] = 1
                data['cn2'] = 2
        elif len(self.sources.keys()) == 2:
            # Assume loss of heterozygosity, since DEL and DUP overlap
            data['status'] = 'neutral'
            data['cn'] = 2
            data['cn1'] = 0
            data['cn2'] = 2

        return pd.Series(data)


def read_cnvkit_no_loh_csv(dict_data, caller):
    """
    Read No-LOH CNVkit CSV file into a Pandas DataFrame.

    :param dict_data: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of No-LOH CNVkit CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'
    path = dict_data[0]
    if dict_data[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)
    df['caller'] = caller

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Set status based on cn
    cond = df['cn'] < ploidy
    df.loc[cond, 'status'] = 'loss'
    cond = df['cn'] == ploidy
    df.loc[cond, 'status'] = 'neutral'
    cond = df['cn'] > ploidy
    df.loc[cond, 'status'] = 'gain'

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


def read_cnvkit_csv(data_dict, caller):
    """
    Read CNVkit CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of CNVkit CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'

    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)
    df['caller'] = caller

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Set status based on cn
    cond = df['cn'] < ploidy
    df.loc[cond, 'status'] = 'loss'
    cond = df['cn'] == ploidy
    df.loc[cond, 'status'] = 'neutral'
    cond = df['cn'] > ploidy
    df.loc[cond, 'status'] = 'gain'

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


def read_purecn_no_loh_csv(data_dict, caller):
    """
    Read No-LOH PureCN CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of No-LOH PureCN CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'

    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)
    df['caller'] = caller

    # Rename columns to match CNV_COLUMNS names
    purecn_no_loh_mapper = {'chrom': 'chromosome', 'loc.start': 'start',
                            'loc.end': 'end', 'C': 'cn'}
    df.rename(mapper=purecn_no_loh_mapper, axis='columns', inplace=True)

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Set status based on cn
    cond = df['cn'] < ploidy
    df.loc[cond, 'status'] = 'loss'
    cond = df['cn'] == ploidy
    df.loc[cond, 'status'] = 'neutral'
    cond = df['cn'] > ploidy
    df.loc[cond, 'status'] = 'gain'

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


def read_purecn_csv(data_dict, caller):
    """
    Read PureCN CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of PureCN CSV file
    """
    skiprows = 0
    header = 0
    sep = ','
    comment = '#'

    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # for now, there is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)
    df['caller'] = caller

    # Rename columns to match CNV_COLUMNS names
    purecn_mapper = {'chr': 'chromosome', 'C': 'cn', 'M': 'cn1'}
    df.rename(mapper=purecn_mapper, axis='columns', inplace=True)

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Set status based on cn
    cond = df['cn'] < ploidy
    df.loc[cond, 'status'] = 'loss'
    cond = df['cn'] == ploidy
    df.loc[cond, 'status'] = 'neutral'
    cond = df['cn'] > ploidy
    df.loc[cond, 'status'] = 'gain'

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    # Calculate cn2
    df['cn2'] = df.cn - df.cn1

    return df, ploidy


def read_sclust(data_dict, caller):
    """
    Read Sclust Allelic States file into a Pandas DataFrame

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :param for_benchmark: Removes excess columns when set to True
    :return: Pandas DataFrame representation of SCLUST Allelic States file
    """

    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'

    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = read_sclust_ploidy(data_dict[1])

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows,
                     header=header, comment=comment, sep=sep)
    df['caller'] = caller

    # Rename columns to match CNV_COLUMNS names
    sclust_mapper = {'Chromosome': 'chromosome', 'CopyNr': 'cn',
                     'A': 'cn1', 'B': 'cn2', 'Start': 'start',
                     'End': 'end'}
    df.rename(mapper=sclust_mapper, axis='columns', inplace=True)

    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Set status based on cn
    cond = df['cn'] < ploidy
    df.loc[cond, 'status'] = 'loss'
    cond = df['cn'] == ploidy
    df.loc[cond, 'status'] = 'neutral'
    cond = df['cn'] > ploidy
    df.loc[cond, 'status'] = 'gain'

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


def read_purple_cnv_somatic(data_dict, caller):
    """
    Read gridss germline csv file

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of Purple CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'

    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = read_purple_ploidy(data_dict[1])

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows,
                     header=header, comment=comment, sep=sep)
    df['caller'] = caller

    # Rename columns to match CNV_COLUMNS names
    purple_mapper = {'chromosome': 'chromosome', 'start': 'start',
                     'end': 'end', 'copyNumber': 'cn',
                     'majorAllelePloidy': 'cn1', 'minorAllelePloidy': 'cn2'}
    df.rename(mapper=purple_mapper, axis='columns', inplace=True)

    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    df = df.round({'cn': 2, 'cn1': 2, 'cn2': 2})

    # Set status based on cn
    cond = df['cn'] < ploidy
    df.loc[cond, 'status'] = 'loss'
    cond = df['cn'] == ploidy
    df.loc[cond, 'status'] = 'neutral'
    cond = df['cn'] > ploidy
    df.loc[cond, 'status'] = 'gain'

    # Set status based on cn
    df['status'] = df.apply(lambda row: set_purple_status(row, ploidy), axis=1)

    # dodati cr kolonu vec ovde

    # Dealing with negative values
    cond = df['cn'] < 0
    df.loc[cond, 'cn'] = 0
    cond = df['cn1'] < 0
    df.loc[cond, 'cn1'] = 0
    cond = df['cn2'] < 0
    df.loc[cond, 'cn2'] = 0

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


def set_purple_additional_status(row, ploidy):
    cn_value = row.loc['cn']
    if cn_value < 0.5:
        status = 'loss_homozygous'
    elif 0.5 <= cn_value < (0.6 * ploidy):
        status = 'loss_heterozygous'
    elif (0.6 * ploidy) <= cn_value < (1.4 * ploidy):
        status = 'neutral'
    elif (1.4 * ploidy) <= cn_value < (3 * ploidy):
        status = 'gain'
    else:
        status = 'gain_high'
    return status


def set_purple_status(row, ploidy):
    cn_value = row.loc['cn']
    if cn_value < 0.5:
        status = 'loss'
    elif 0.5 <= cn_value < (0.6 * ploidy):
        status = 'loss'
    elif (0.6 * ploidy) <= cn_value < (1.4 * ploidy):
        status = 'neutral'
    elif (1.4 * ploidy) <= cn_value < (3 * ploidy):
        status = 'gain'
    else:
        status = 'gain'
    return status


def read_controlfreec_csv(data_dict, caller, for_benchmark=True):
    """
    Read ControlFREEC CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :param for_benchmark: Removes excess columns when set to True
    :return: Pandas DataFrame representation of ControlFREEC CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'

    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = read_controlfreec_ploidy(data_dict[1])

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)
    df['caller'] = caller

    # Rename columns to match CNV_COLUMNS names
    control_freec_mapper = {'chr': 'chromosome', 'copy number': 'cn'}
    df.rename(mapper=control_freec_mapper, axis='columns', inplace=True)

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Fill genotype with unknown if missing
    if 'genotype' not in df.columns:
        df['genotype'] = '-'

    # Infer cn1 and cn2 from genotype field
    df[['cn1', 'cn2']] = df['genotype'].apply(extract_cn1_and_cn2)

    # df['status'] = df.apply(
    # lambda row: set_controlfreec_status(row, ploidy), axis=1
    # )

    # Set DTypes
    if for_benchmark:
        df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


def set_controlfreec_status(row, ploidy):
    cn1_value = row.loc['cn1']
    cn2_value = row.loc['cn2']
    freec_status = row.loc['status']
    if freec_status == 'loss':
        if cn1_value == cn2_value:
            status = 'loss_homozygous'
        elif cn1_value != cn2_value:
            status = 'loss_heterozygous'
        else:
            status = freec_status
    else:
        status = freec_status

    return status


def read_controlfreec_no_loh_csv(data_dict, caller):
    """
    Read ControlFREEC no-LOH CSV file into a Pandas DataFrame.

    :param data_dict: Path to ControlFREEC no-LOH CSV file
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of ControlFREEC no-LOH CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'

    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)
    df['caller'] = caller

    # Rename columns to match CNV_COLUMNS names
    control_freec_mapper = {'chr': 'chromosome', 'copy number': 'cn'}
    df.rename(mapper=control_freec_mapper, axis='columns', inplace=True)

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Fill genotype with unknown if missing
    if 'genotype' not in df.columns:
        df['genotype'] = '-'

    # Infer cn1 and cn2 from genotype field
    df[['cn1', 'cn2']] = df['genotype'].apply(extract_cn1_and_cn2)

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df


def read_controlfreec_no_header_csv(data_dict, caller, for_benchmark=True):
    """
    Read ControlFREEC no-header CSV file into a Pandas DataFrame.

    :param data_dict: Path to ControlFREEC no-header CSV file
    :param caller: Name of the caller
    :param for_benchmark: Removes excess columns when set to True
    :return: Pandas DataFrame representation of ControlFREEC no-header CSV file
    """
    skiprows = 0
    header = None
    sep = '\t'
    comment = '#'

    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)

    df.columns = ['chromosome', 'start', 'end', 'cn',
                  'status', 'genotype', 'uncertainty']
    df['caller'] = caller

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Fill genotype with unknown if missing
    if 'genotype' not in df.columns:
        df['genotype'] = '-'

    # Infer cn1 and cn2 from genotype field
    df[['cn1', 'cn2']] = df['genotype'].apply(extract_cn1_and_cn2)

    # Set DTypes
    if for_benchmark:
        df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df


def read_controlfreec_no_loh_no_header_csv(data_dict, caller):
    """
    Read ControlFREEC no-LOH no-header CSV file into a Pandas DataFrame.

    :param data_dict: Path to ControlFREEC CSV file
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of ControlFREEC
             no-LOH no-header CSV file
    """
    skiprows = 0
    header = None
    sep = '\t'
    comment = '#'

    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)

    df.columns = ['chromosome', 'start', 'end', 'cn', 'status']
    df['caller'] = caller

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Fill genotype with unknown if missing
    if 'genotype' not in df.columns:
        df['genotype'] = '-'

    # Infer cn1 and cn2 from genotype field
    df[['cn1', 'cn2']] = df['genotype'].apply(extract_cn1_and_cn2)

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df


def read_cnvnator_csv(data_dict, caller):
    """
    Read CNVnator CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of CNVnator CSV file
    """
    skiprows = 0
    header = None
    sep = '\t'
    comment = '#'
    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    cnvnator_columns = ['call', 'coordinates', 'CNV_size', 'normalized_RD',
                        'e-val1', 'e-val2', 'e-val3', 'e-val4', 'q0']
    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, names=cnvnator_columns, comment=comment)
    df['caller'] = caller

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Convert CNVnator columns to CNV columns
    df['chromosome'] = df['coordinates'].apply(
        lambda x: ':'.join(x.split(':')[:-1]))
    df['start'] = df['coordinates'].apply(
        lambda x: (x.split(':')[-1]).split('-')[0])
    df['end'] = df['coordinates'].apply(
        lambda x: (x.split(':')[-1]).split('-')[-1])
    df['status'] = df['call'].map({'deletion': 'loss', 'duplication': 'gain'})

    #Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


def read_purple_tsv(path, caller):
    """
    Read PURPLE TSV file into a Pandas DataFrame.

    :param path: Path to PURPLE TSV file
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of PURPLE TSV file
    """
    skiprows = 0
    header = None
    sep = '\t'
    comment = '#'

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)
    df = df.iloc[1:]
    df.columns = ['chromosome', 'start', 'end', 'copyNumber',
                  'bafCount', 'observedBAF', 'baf', 'segmentStartSupport',
                  'segmentEndSupport', 'method', 'depthWindowCount',
                  'gcContent',
                  'minStart', 'maxStart', 'minorAllelePloidy',
                  'majorAllelePloidy']
    df['caller'] = caller

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Convert PURPLE columns to CNV columns
    df["copyNumber"] = pd.to_numeric(df["copyNumber"], downcast="float")
    df['status'] = ['gain' if x > 2 else 'loss' for x in df['copyNumber']]
    df.loc[df['copyNumber'] == 2, 'status'] = 'neutral'

    return df


def read_icr96_truth_csv(data_dict, caller):
    """
    Read ICR96 truth CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of ICR96 truth CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'
    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)
    df['caller'] = caller

    # Rename columns to match CNV_COLUMNS names
    icr96_mapper = {'ExonCNVType': 'cn'}
    df.rename(mapper=icr96_mapper, axis='columns', inplace=True)

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Set status based on cn
    cond = df['cn'] < ploidy
    df.loc[cond, 'status'] = 'loss'
    cond = df['cn'] == ploidy
    df.loc[cond, 'status'] = 'neutral'
    cond = df['cn'] > ploidy
    df.loc[cond, 'status'] = 'gain'

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


def read_sequenza_csv(data_dict, caller):
    """
    Read Sequenza CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of Sequenza CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'
    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)
    df['caller'] = caller

    # Rename columns to match CNV_COLUMNS names
    sequenza_mapper = {'start.pos': 'start', 'end.pos': 'end', 'CNt': 'cn',
                       'A': 'cn1', 'B': 'cn2'}
    df.rename(mapper=sequenza_mapper, axis='columns', inplace=True)

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Set status based on cn
    cond = df['cn'] < ploidy
    df.loc[cond, 'status'] = 'loss'
    cond = df['cn'] == ploidy
    df.loc[cond, 'status'] = 'neutral'
    cond = df['cn'] > ploidy
    df.loc[cond, 'status'] = 'gain'

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


def read_gatk_called_csv(data_dict, caller, for_benchmark=True):
    """
    Read GATK4 CNV CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :param for_benchmark: Removes excess columns when set to True
    :return: Pandas DataFrame representation of GATK4 CNV CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '@'

    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)
    df['caller'] = caller

    # Rename columns to match CNV_COLUMNS names
    gatk_called_mapper = {'CONTIG': 'chromosome', 'START': 'start',
                          'END': 'end'}
    df.rename(mapper=gatk_called_mapper, axis='columns', inplace=True)

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Set status based on CALL column: '-': 'loss', '0': 'neutral', '+': 'gain'
    df['status'] = df['CALL'].map({'-': 'loss', '0': 'neutral', '+': 'gain'})

    # Set cn based on MEAN_LOG2_COPY_RATIO column, assume 2 is baseline
    df['cn'] = (2 ** (df['MEAN_LOG2_COPY_RATIO'] + 1)).round(1)

    # Set DTypes
    if for_benchmark:
        df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


def read_precise_and_majority_calls(data_dict, caller, for_benchmark=True):
    """
    Read precise / majority file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :param for_benchmark: Removes excess columns when set to True
    :return: Pandas DataFrame representation of GATK4 CNV CSV file
    """
    #     skiprows = 0
    #     header = 0
    #     sep = '\t'
    #     comment = '@'

    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path)
    df['caller'] = caller

    # Rename columns to match CNV_COLUMNS names
    majority_called_mapper = {'chromosome': 'chromosome', 'start': 'start',
                              'end': 'end', 'final_call': 'status'}
    df.rename(mapper=majority_called_mapper, axis='columns', inplace=True)

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Set DTypes
    if for_benchmark:
        df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


def read_gatk_model_csv(data_dict, caller, for_benchmark=True):
    """
    Read GATK Model Final SEG CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :param for_benchmark: Removes excess columns when set to True
    :return: Pandas DataFrame representation of GATK Model Final SEG CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '@'
    na_values = 'NaN'

    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment, na_values=na_values)
    df['caller'] = caller

    # Rename columns to match CNV_COLUMNS names
    gatk_model_mapper = {'CONTIG': 'chromosome', 'START': 'start',
                         'END': 'end'}
    df.rename(mapper=gatk_model_mapper, axis='columns', inplace=True)

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Set DTypes
    if for_benchmark:
        df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


def read_facets_csv(data_dict, caller):
    """
    Read Facets CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of Facets CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'
    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = read_facets_ploidy(data_dict[1])

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)
    df['caller'] = caller

    # Rename columns to match CNV_COLUMNS names
    facets_cnv_mapper = {'chrom': 'chromosome', 'tcn.em': 'cn',
                         'lcn.em': 'cn1'}
    df.rename(mapper=facets_cnv_mapper, axis='columns', inplace=True)

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Set status based on cn
    cond = df['cn'] < ploidy
    df.loc[cond, 'status'] = 'loss'
    cond = df['cn'] == ploidy
    df.loc[cond, 'status'] = 'neutral'
    cond = df['cn'] > ploidy
    df.loc[cond, 'status'] = 'gain'
    # df['status'] = df.apply(lambda row: set_range_status(row, ploidy), axis=1)

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    # Set cn2 based on cn1 column
    df['cn2'] = df['cn'] - df['cn1']

    return df, ploidy


def read_scnvsim_csv(data_dict, caller):
    """
    Read SCNVSim simulated CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of SCNVSim simulated CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'
    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)
    df['caller'] = caller

    # Rename columns to match CNV_COLUMNS names
    scnvsim_mapper = {'chr': 'chromosome', 'Start': 'start', 'End': 'end',
                      'Copy_Number': 'cn'}
    df.rename(mapper=scnvsim_mapper, axis='columns', inplace=True)

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    # SCNVSim reports difference from baseline (2)
    df['cn'] = df['cn'] + 2

    # Set status based on cn
    cond = df['cn'] < 2
    df.loc[cond, 'status'] = 'loss'
    cond = df['cn'] == 2
    df.loc[cond, 'status'] = 'neutral'
    cond = df['cn'] > 2
    df.loc[cond, 'status'] = 'gain'

    return df, ploidy


def read_simulatecnvs_csv(data_dict, caller):
    """
    Read SimulateCNVs simulated CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of SimulateCNVs simulated CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'
    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)
    df['caller'] = caller

    # Rename columns to match CNV_COLUMNS names
    simulatecnvs_mapper = {'chr': 'chromosome', 'copy_number': 'cn'}
    df.rename(mapper=simulatecnvs_mapper, axis='columns', inplace=True)

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Set status based on cn
    cond = df['cn'] < ploidy
    df.loc[cond, 'status'] = 'loss'
    cond = df['cn'] == ploidy
    df.loc[cond, 'status'] = 'neutral'
    cond = df['cn'] > ploidy
    df.loc[cond, 'status'] = 'gain'

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


def read_varsimlab_csv(data_dict, caller):
    """
    Read VarSimLab simulated CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of VarSimLab simulated CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'
    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment)
    df['caller'] = caller

    # Rename columns to match CNV_COLUMNS names
    varsimlab_mapper = {'chr': 'chromosome', 'location': 'start',
                        'copies1': 'cn1', 'copies2': 'cn2'}
    df.rename(mapper=varsimlab_mapper, axis='columns', inplace=True)

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Strip chromosome column
    df['chromosome'] = df['chromosome'].str.strip()

    # Set end based on start and seq_size columns
    df['end'] = df['start'] + df['seq_size']

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    # Set cn based on cn1 and cn2 columns
    df['cn'] = df['cn1'] + df['cn2']

    # Set status based on cn
    cond = df['cn'] < ploidy
    df.loc[cond, 'status'] = 'loss'
    cond = df['cn'] == ploidy
    df.loc[cond, 'status'] = 'neutral'
    cond = df['cn'] > ploidy
    df.loc[cond, 'status'] = 'gain'

    return df, ploidy


def read_bed_csv(data_dict, caller):
    """
    Read BED CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of BED CSV file
    """
    skiprows = 0
    header = None
    sep = '\t'
    comment = '#'
    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    bed_columns = ['chromosome', 'start', 'end']
    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, names=bed_columns, comment=comment)
    df['caller'] = caller

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


def read_sbg_conseca_cnv(data_dict, caller, for_benchmark=True):
    """
    Read SBG Conseca CNV CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :param for_benchmark: Removes excess columns when set to True
    :return: Pandas DataFrame representation of SBG Conseca CNV CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'
    na_values = '.'
    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     comment=comment, na_values=na_values, sep=sep)
    df['caller'] = caller

    if for_benchmark:
        # Rename columns to match CNV_COLUMNS names
        conseca_mapper = {'freec_cn': 'cn', 'freec_cn1': 'cn1',
                          'freec_cn2': 'cn2'}
        df.rename(mapper=conseca_mapper, axis='columns', inplace=True)

    # Add missing columns and set them to np.NaN
    for column in CNV_COLUMNS:
        if column not in df.columns:
            df[column] = np.NaN

    if for_benchmark:
        # Filter out low and no confidence calls
        cond = df['confidence'] == 'high'
        df = df[cond]

    # Set DTypes
    if for_benchmark:
        df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


def read_sveto_prep_cnv_csv(data_dict, caller):
    """
    Read SVeto Prepare CNV CSV file into a Pandas DataFrame.

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :return: Pandas DataFrame representation of SVeto Prepare CNV CSV file
    """
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'
    path = data_dict[0]
    if data_dict[1] is None:
        ploidy = 2
    else:
        ploidy = 2  # There is no calculation of ploidy for this caller

    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows, header=header,
                     sep=sep, comment=comment, na_values='.')
    df['caller'] = caller

    # Set DTypes
    df = df[CNV_COLUMNS].copy()
    for key in CNV_DTYPES:
        df[key] = df[key].astype(CNV_DTYPES[key])

    return df, ploidy


# TODO sync caller names, formats and file extensions
def read_csv(data_dict, caller, cnv_format=None):
    """
    Read CNV CSV with the appropriate reader, based on the file extension.

    *.no_loh.cns:       read_cnvkit_no_loh_csv()
    *.cns:              read_cnvkit_csv()
    *.called.seg:       read_gatk_called_csv()
    *.modelFinal.seg:   read_gatk_model_csv()
    *.seg:              read_purecn_no_loh_csv()
    *_loh.csv:          read_purecn_csv()
    *.value.txt:        read_controlfreec_csv()
    *.CNV_results.txt:  read_varsimlab_csv()
    *_results.txt:      read_cnvnator_csv()
    *.exon.csv:         read_icr96_truth_csv()
    *segments.txt:      read_sequenza_csv()
    *.cncf.tsv:         read_facets_csv()
    *.scnvsim.bed:      read_scnvsim_csv()
    *_exon.bed:         read_simulatecnvs_csv()
    *.bed:              read_bed_csv()
    *.sveto-prep.cnv:   read_sveto_prep_cnv_csv()
    *.conseca.tsv:      read_sbg_conseca_cnv()
    *.cnv.somatic.tsv:  read_purple_cnv_somatic()

    :param data_dict: Dict with path to data file and ploidy info
    :param caller: Name of the caller
    :param cnv_format: Format of the CNV file
    :return: Pandas DataFrame representation of CNV CSV file
    """

    cnv_readers = {'bed': read_bed_csv,
                   'cnvkit': read_cnvkit_csv,
                   'cnvkit_no_loh': read_cnvkit_no_loh_csv,
                   'cnvnator': read_cnvnator_csv,
                   'controlfreec': read_controlfreec_csv,
                   'controlfreec_no_loh': read_controlfreec_no_loh_csv,
                   'controlfreec_no_header': read_controlfreec_no_header_csv,
                   'controlfreec_no_loh_no_header':
                       read_controlfreec_no_loh_no_header_csv,
                   'facets': read_facets_csv,
                   'gatk_called': read_gatk_called_csv,
                   'gatk_model': read_gatk_model_csv,
                   'icr96': read_icr96_truth_csv,
                   'purecn': read_purecn_csv,
                   'purecn_no_loh': read_purecn_no_loh_csv,
                   'scnvsim': read_scnvsim_csv,
                   'sequenza': read_sequenza_csv,
                   'simulatecnvs': read_simulatecnvs_csv,
                   'varsimlab': read_varsimlab_csv,
                   'sveto_prep_cnv': read_sveto_prep_cnv_csv,
                   'sbg_conseca_cnv': read_sbg_conseca_cnv,
                   'sclust': read_sclust,
                   'purple_somatic': read_purple_cnv_somatic,
                   'majority_calls': read_precise_and_majority_calls,
                   'precise_calls': read_precise_and_majority_calls
                   }
    path = data_dict[0]

    if cnv_format in cnv_readers:
        df = cnv_readers[cnv_format](path, caller)
        df = df.sort_values(['chromosome', 'start'])
        return df
    # CNVkit no LOH
    if path.endswith('.no_loh.cns') or path.endswith('_cnv.txt'):
        df, ploidy = read_cnvkit_no_loh_csv(data_dict, caller)
    # CNVkit
    elif path.endswith('.cns'):
        df, ploidy = read_cnvkit_csv(data_dict, caller)
    # GATK called
    elif path.endswith('.called.seg'):
        df, ploidy = read_gatk_called_csv(data_dict, caller)
    # GATK model
    elif path.endswith('.modelFinal.seg'):
        df, ploidy = read_gatk_model_csv(data_dict, caller)
    # PureCN no LOH
    elif path.endswith('.seg'):
        df, ploidy = read_purecn_no_loh_csv(data_dict, caller)
    # PureCN
    elif path.endswith('_loh.csv'):
        df, ploidy = read_purecn_csv(data_dict, caller)
    # ControlFREEC
    elif path.endswith('.value.txt'):
        df, ploidy = read_controlfreec_csv(data_dict, caller)
    # VarSimLab
    elif path.endswith('.CNV_results.txt'):
        df, ploidy = read_varsimlab_csv(data_dict, caller)
    # CNVnator
    elif path.endswith('calling_results.txt'):
        df, ploidy = read_cnvnator_csv(data_dict, caller)
    # ICR96 truth
    elif path.endswith('.exon.csv'):
        df, ploidy = read_icr96_truth_csv(data_dict, caller)
    # Sequenza
    elif path.endswith('.segments.txt'):
        df, ploidy = read_sequenza_csv(data_dict, caller)
    # Facets
    elif path.endswith('.cncf.tsv'):
        df, ploidy = read_facets_csv(data_dict, caller)
    # SCNVSim
    elif path.endswith('.scnvsim.bed'):
        df, ploidy = read_scnvsim_csv(data_dict, caller)
    # SimulateCNVs
    elif path.endswith('.overlap_exon.bed'):
        df, ploidy = read_simulatecnvs_csv(data_dict, caller)
    # BED
    elif path.endswith('.bed'):
        df, ploidy = read_bed_csv(data_dict, caller)
    # SVeto Prepare SV -> CNV
    elif path.endswith('.sveto-prep.cnv'):
        df, ploidy = read_sveto_prep_cnv_csv(data_dict, caller)
    # SBG Conseca CNV
    elif path.endswith('.conseca.tsv'):
        df, ploidy = read_sbg_conseca_cnv(data_dict, caller)
    # Sclust
    elif path.endswith('_allelic_states.txt'):
        df, ploidy = read_sclust(data_dict, caller)
    # Purple Somatic
    elif path.endswith('.cnv.somatic.tsv'):
        df, ploidy = read_purple_cnv_somatic(data_dict, caller)
    elif path.endswith('.majority_calls.csv'):
        df, ploidy = read_precise_and_majority_calls(data_dict, caller)
    elif path.endswith('.precise_calls.csv'):
        df, ploidy = read_precise_and_majority_calls(data_dict, caller)
    else:
        extension = 'unknown'
        if path:
            extension = path.split('.')[-1]
        raise UserWarning(
            'File type "{} {}" is not supported.'.format(extension, path))

    df = df.sort_values(['chromosome', 'start'])
    return df, ploidy


def read_ploidy(path, cnv_format=None):
    """
    Read CNV CSV with the appropriate reader, based on the file extension.

    *.out_cn_summary.txt: read_sclust_ploidy()
    *.CNV-loh.config.txt: read_controlfreec_ploidy()
    *.purple.purity.tsv:  read_purple_ploidy()

    :param path: Path to CNV ploidy file
    :param cnv_format: Format of the ploidy file
    :return: Pandas DataFrame representation of CNV CSV file
    """

    ploidy_readers = {
        'sclust': read_sclust_ploidy,
        'controlfreec': read_controlfreec_ploidy,
        'purple': read_purple_ploidy
    }

    # SCLUST
    if path.endswith('.out_cn_summary.txt'):
        ploidy = read_sclust_ploidy(path)
    # Control FREEC
    elif path.endswith('.CNV-loh.config.txt'):
        ploidy = read_controlfreec_ploidy(path)
    # PURPLE
    elif path.endswith('.purple.purity.tsv'):
        ploidy = read_purple_ploidy(path)
    else:
        extension = 'unknown'
        if path:
            extension = path.split('.')[-1]
        raise UserWarning(
            'File type "{} {}" is not supported.'.format(extension, path))

    return ploidy


def read_purple_ploidy(path):
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'
    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows,
                     header=header, comment=comment, sep=sep,
                     dtype={"ploidy": float, "purity": float})
    ploidy = df['ploidy'][0]
    ploidy_to_return = round(float(ploidy), 2)

    return ploidy_to_return


def read_controlfreec_ploidy(path):
    f = open(path, "r")
    ploidy = 2
    with open(path) as openfile:
        for line in openfile:
            if 'Output_Ploidy' in line:
                parts = line.split('\t')
                ploidy = float(parts[1].split('\n')[0])
                return ploidy
    return ploidy


def read_sclust_ploidy(path):
    f = open(path, "r")
    header_line = f.readline()
    header_indexes = header_line.split('\t')
    value_line = f.readline()
    value_indexes = value_line.split('\t')
    data_dict = dict(zip(header_indexes, value_indexes))
    ploidy_value = (data_dict['ploidy'])
    ploidy = float(ploidy_value)
    ploidy = int(ploidy)
    return ploidy


def read_facets_ploidy(path):
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'
    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows,
                     header=header, comment=comment, sep=sep)
    ploidy = df['ploidy'][0]
    ploidy_to_return = round(float(ploidy), 2)

    return ploidy_to_return


def convert_cnv_to_ss(df):
    """
    Output CN calls from CN Dataframe to ShatterSeek SSCN calls.

    :param df: CN DataFrame
    :return: Dataframe in format of ShatterSeek CN (SSCN)
    """
    ss_mapper = {'chromosome': 'chrom', 'cn': 'total_cn'}
    df_ss = df[['chromosome', 'start', 'end', 'cn']].copy()
    df_ss.rename(mapper=ss_mapper, axis='columns', inplace=True)

    # Remove CHR, Chr, chr from chromosome name
    # This also ensures that X is X and not x
    df_ss.chrom = df_ss.chrom.str.upper().str.replace('CHR', '')

    # Keep 1, 2, ... 22, X chromosomes
    valid_chroms = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X'
    valid_chroms = valid_chroms.split()
    cond = df_ss.chrom.isin(valid_chroms)
    df_ss = df_ss[cond]

    # Set dtypes for pos1 and pos2
    df_ss.start = df_ss.start.astype(int)
    df_ss.end = df_ss.end.astype(int)

    return df_ss


def write_cnv(df, columns, output):
    """
    Write CN DataFrame to provided output stream.

    :param df: DataFrame to output
    :param columns: List of columns to output (also sets column order)
    :param output: Stream (opened file) to write to
    :return:
    """
    cnv = df[columns].to_csv(sep='\t', index=False, na_rep='.')
    output.write(cnv)


def fix_chromosomes(df):
    """
    Remove CHR/Chr/chr from contig ID and filter out alternative contigs - leave
    only 1, 2, ..., 22, X and Y contigs.

    :param df: DataFrame to be processed, must contain 'chromosome' column
    :return: Processed and filtered DataFrame
    """
    # Remove CHR, Chr, chr from chromosome name
    # This also ensures that X is X and not x
    df['chromosome'] = df['chromosome'].str.upper().str.replace('CHR', '')

    # Keep 1, 2, ... 22, X, Y chromosomes
    chromosomes = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y'
    chromosomes = chromosomes.split()
    cond = df['chromosome'].isin(chromosomes)

    return df[cond]


def extract_cn1_and_cn2(genotype):
    """
    Infer cn1 and cn2 values from Control-FREEC 'genotype' column
    :param genotype: Control-FREEC 'genotype' column (e.g. 'AAABB')
    :return: Series object with inferred 'cn1' and 'cn2' columns
    """
    cnt_a = genotype.count('A')
    cnt_b = genotype.count('B')

    if cnt_a == 0 and cnt_b == 0:
        return pd.Series({'cn1': np.NaN, 'cn2': np.NaN})
    else:
        return pd.Series({'cn1': cnt_a, 'cn2': cnt_b})


def remove_overlaps_from_df(df, method='SWAP'):
    """
    Look for overlaps in a DataFrame (region[i+1].start < region[i].end) and
    fix it by swapping those 2 values.

    :param df: DataFrame to remove overlaps from
    :param method: SWAP [n+1].start and [n].end; SPLIT set both to rounded mean
    :return: DataFrame with no overlaps
    """
    cur = None
    for index, row in df.sort_values(
            ['chromosome', 'start', 'end']).iterrows():
        if cur is None:
            cur = {'index': index, 'chromosome': row.chromosome,
                   'start': row.start, 'end': row.end}
            continue

        if cur['chromosome'] == row.chromosome and cur['end'] > row.start:
            if cur['end'] >= row.end:
                raise UserWarning(
                    'Error removing overlap for {} and {}'.format(cur, row))

            left = row.start
            right = cur['end']
            if method == 'SWAP':
                print(
                    'Warning: Removing self-intersection at {}({}:{})'.format(
                        row.chromosome, int(left), int(right)))

            if method == 'SWAP':
                df.loc[cur['index'], 'end'] = left
                df.loc[index, 'start'] = right
            elif method == 'SPLIT':
                middle = int((left + right) / 2)
                df.loc[cur['index'], 'end'] = middle
                df.loc[index, 'start'] = middle

            cur = {'index': index, 'chromosome': row.chromosome,
                   'start': right, 'end': row.end}
        else:
            cur = {'index': index, 'chromosome': row.chromosome,
                   'start': row.start, 'end': row.end}


def to_regions(dataframe):
    """
    Create a list of Regions based on a DataFrame.

    :param dataframe: DataFrame to create Regions list from
    :return: List of Regions corresponding to the input DataFrame
    """
    return [Region(r.chromosome, r.start, r.end, [{r.caller: r}])
            for _, r in dataframe.iterrows()]


def combine_regions(one_regions, two_regions, allow_duplicates=False):
    """
    Process 'one' and 'two' Region lists to generate combined Region list.

    :param one_regions: First list of regions
    :param two_regions: Second list of regions
    :param allow_duplicates: Duplicate callers allowed in the region
    :return: List of combined regions
    """
    tail = None
    regions = []
    for region in sorted(one_regions + two_regions, key=repr):
        if tail is None:
            tail = region
        else:
            outs, tail = tail.intersect(region, allow_duplicates)
            regions.extend(outs)

    if tail is not None:
        regions.append(tail)

    return regions


def filter_regions(regions, caller):
    """
    Filter out Region objects which don't contain caller in sources. Allows use
    of a BED file for filtering.

    :param regions: List of Regions to filter
    :param caller: Caller needed for Region to be in output list
    :return: List of Regions with caller present in sources
    """
    return [r for r in regions if caller in r.sources]


def to_dataframe(regions, callers):
    """
    Create a DataFrame from Regions list for given analysis and callers.

    :param regions: List of Regions to process
    :param callers: List of callers to take into account
    :return: DataFrame with Region coordinates and data for each caller
    """
    return pd.DataFrame((r.to_series(callers) for r in regions))


def to_dataframe_with_cn(regions, callers):
    """
    Create a DataFrame from Regions list for given analysis and callers.

    :param regions: List of Regions to process
    :param callers: List of callers to take into account
    :return: DataFrame with Region coordinates and data for each caller
    """
    return pd.DataFrame((r.to_series_with_cn(callers) for r in regions))


def to_dataframe_with_cr(regions, callers_and_ploidy_dict):
    """
    Create a DataFrame from Regions list for given analysis and callers.

    :param regions: List of Regions to process
    :param callers_and_ploidy_dict: List of callers and ploidy info
    :return: DataFrame with Region coordinates and data for each caller
    """
    return pd.DataFrame(
        (r.to_series_with_cr(callers_and_ploidy_dict) for r in regions))


def to_dataframe_for_consensus(regions):
    """
    Create a DataFrame from Regions list for consensus calling.

    :param regions: List of Regions to process
    :return: DataFrame with Region coordinates and consensus info
    """
    columns = ['chromosome', 'start', 'end', 'status', 'confidence', 'support',
               'freec_status', 'freec_genotype', 'freec_cn', 'freec_cn1',
               'freec_cn2', 'gatk_status', 'gatk_cn', 'gatk_call',
               'gatk_mean_log2_copy_ratio',
               'gatk_minor_allele_fraction_posterior_50']
    return pd.DataFrame((r.to_series_for_consensus() for r in regions))[
        columns]


def to_dataframe_for_gatk(regions):
    """
    Create a DataFrame from Regions list for GATK cnLOH.

    :param regions: List of Regions to process
    :return: DataFrame with Region coordinates and GATK cnLOH info
    """
    columns = ['chromosome', 'start', 'end', 'length', 'status', 'cn', 'cn1',
               'cn2', 'call', 'mean_log2_copy_ratio',
               'log2_copy_ratio_posterior_50',
               'minor_allele_fraction_posterior50']
    return pd.DataFrame((r.to_series_for_gatk() for r in regions))[columns]


def to_dataframe_for_sv(regions):
    """
    Create a DataFrame from Regions list initiated from SVs.

    :param regions: List of Regions to process
    :return: DataFrame with Region coordinates and data for each caller
    """
    return pd.DataFrame((r.to_series_for_sv() for r in regions))


def calculate_metrics(regions, tru_caller, que_caller, custom_metrics=False):
    """
    Calculate metrics for given DataFrame.

    :param regions: Regions to evaluate
    :param tru_caller: Name of truth caller
    :param que_caller: Name of caller to compare to truth caller
    :param custom_metrics: If selected calculate TP, OE and UE only
    :return: Series containing metrics
    """
    metrics_columns = ['caller', 'total', 'tp', 'tn', 'fp', 'fn', 'oe', 'ue',
                       'recall', 'fscore', 'precision']

    custom_metrics_columns = ['caller', 'total', 'tp', 'oe', 'ue',
                              'tp_percent',
                              'oe_percent', 'ue_percent']

    tru_status = '{}_status'.format(tru_caller)
    que_status = '{}_status'.format(que_caller)

    if custom_metrics:
        # Replace NaN in truth status with neutral
        cond = regions[tru_status].isnull()
        regions.loc[cond, tru_status] = 'neutral'

        # Replace NaN in query status with neutral
        cond = regions[que_status].isnull()
        regions.loc[cond, que_status] = 'neutral'

        cond_tp = regions[tru_status] == regions[que_status]
        tp = regions[cond_tp].length.sum()

        cond_tru_loss = regions[tru_status].astype(str) == 'loss'
        cond_tru_neut = regions[tru_status].astype(str) == 'neutral'
        cond_tru_gain = regions[tru_status].astype(str) == 'gain'
        cond_que_loss = regions[que_status].astype(str) == 'loss'
        cond_que_neut = regions[que_status].astype(str) == 'neutral'
        cond_que_gain = regions[que_status].astype(str) == 'gain'

        cond_oe = (cond_tru_loss & cond_que_neut) | (
                cond_tru_loss & cond_que_gain) | (
                              cond_tru_neut & cond_que_gain)
        oe = regions[cond_oe].length.sum()

        cond_ue = (cond_tru_gain & cond_que_neut) | (
                cond_tru_gain & cond_que_loss) | (
                              cond_tru_neut & cond_que_loss)
        ue = regions[cond_ue].length.sum()

        tp_percent = tp / (tp + oe + ue)
        oe_percent = oe / (tp + oe + ue)
        ue_percent = ue / (tp + oe + ue)

        total = tp + oe + ue

        return pd.Series(
            [que_caller, total, tp, oe, ue, tp_percent, oe_percent,
             ue_percent],
            index=custom_metrics_columns)

    else:
        # Replace neutral in truth status with NaN
        cond = regions[tru_status] == 'neutral'
        regions.loc[cond, tru_status] = np.NaN

        # Replace neutral in query status with NaN
        cond = regions[que_status] == 'neutral'
        regions.loc[cond, que_status] = np.NaN

        cond_tp = regions[tru_status] == regions[que_status]
        tp = regions[cond_tp].length.sum()

        cond_tn = regions[tru_status].isnull() & regions[que_status].isnull()
        tn = regions[cond_tn].length.sum()

        cond_fp = regions[tru_status].isnull() & regions[que_status].notnull()
        fp = regions[cond_fp].length.sum()

        cond_fn = regions[tru_status].notnull() & regions[que_status].isnull()
        fn = regions[cond_fn].length.sum()

        cond_oe = (regions[tru_status].astype(str) == 'loss') & (
                regions[que_status].astype(str) == 'gain')
        oe = regions[cond_oe].length.sum()

        cond_ue = (regions[tru_status].astype(str) == 'gain') & (
                regions[que_status].astype(str) == 'loss')
        ue = regions[cond_ue].length.sum()

        recall = tp / (tp + fn + oe + ue)
        precision = tp / (tp + fp + oe + ue)
        fscore = 0
        if recall + precision != 0:
            fscore = 2 * recall * precision / (recall + precision)

        total = tp + tn + fp + fn + oe + ue

        return pd.Series(
            [que_caller, total, tp, tn, fp, fn, oe, ue, recall, fscore,
             precision], index=metrics_columns)


def calculate_metrics_zare(truth, query, overlap_threshold=0.8):
    """
    Calculate metrics using Zare et al. method for given DataFrame.

    :param truth: Truth DataFrame
    :param query: Query DataFrame
    :param overlap_threshold: Same call if overlap greater than threshold
    :return: Series containing metrics
    """

    def compare_rows(row1, row2, overlap_threshold=0.8):
        """
        Inline helper method.

        :param row1: First row
        :param row2: Second row
        :param overlap_threshold: Same call if overlap greater than threshold
        :return: True if same call and overlap > threshold
        """
        chromosome = row1['chromosome'] == row2['chromosome']
        status = row1['status'] == row2['status']
        overlap = row1['start'] < row2['end'] and row1['end'] > row2['start']
        intersection = min(row1['end'], row2['end']) - max(row1['start'],
                                                           row2['start'])
        union = max(row1['end'], row2['end']) - min(row1['start'],
                                                    row2['start'])
        overlap_enough = intersection / union > overlap_threshold

        if chromosome and status and overlap and overlap_enough:
            return True
        else:
            return False

    zare_metrics_columns = ['caller', 'overlap_threshold', 'tt', 'fn', 'tp',
                            'fp', 'tq', 'recall', 'fscore', 'precision']
    que_caller = query['caller'].iloc[-1]

    # tt - total truth calls
    # tq - total query calls
    # tp - true positives
    # fn - false negatives (fn = tt - tp)
    # fp - false positives (fp = tq - tp)

    # Keep only loss and gain regions
    # Truth
    cond = truth['status'].isin(['loss', 'gain'])
    truth = truth[cond]
    # Query
    cond = query['status'].isin(['loss', 'gain'])
    query = query[cond]

    tt = len(truth)
    tq = len(query)

    tp = sum((compare_rows(truth_row, query_row, overlap_threshold)
              for _, truth_row
              in truth.iterrows()
              for _, query_row
              in query[(query['chromosome'] == truth_row['chromosome']) &
                       (query['status'] == truth_row['status']) &
                       (query['start'] < truth_row['end']) &
                       (query['end'] > truth_row['start'])].iterrows()))

    fn = tt - tp
    fp = tq - tp

    recall = tp / tt
    precision = 0
    if tq != 0:
        precision = tp / tq
    fscore = 0
    if recall + precision != 0:
        fscore = 2 * recall * precision / (recall + precision)

    return pd.Series(
        [que_caller, overlap_threshold, tt, fn, tp, fp, tq, recall, fscore,
         precision], index=zare_metrics_columns)


def powerset(iterable):
    """
    Create all combinations of iterable's elements. Used for Venn analysis.

    Example:
    powerset([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)

    :param iterable: List (iterable) of elements to combine
    :return: Iterator of element combinations tuples
    """
    s = list(iterable)  # allows duplicate elements
    return chain.from_iterable(
        combinations(s, r) for r in range(1, len(s) + 1))


def venn(df, callers):
    """
    Analyze congruence of all caller combinations. Find length of loss, neutral
    and gain regions where each combination 'agrees'.

    :param df: Combined regions DataFrame
    :param callers: List of callers to analyze
    :return: DataFrame containing results
    """
    data = []

    loss_conds = {}
    neut_conds = {}
    gain_conds = {}

    for caller in callers:
        status = '{}_status'.format(caller)
        loss_conds[caller] = df[status].isin(['loss'])
        neut_conds[caller] = df[status].isnull() | df[status].isin(['neutral'])
        gain_conds[caller] = df[status].isin(['gain'])

    for combo in powerset(callers):
        row = {}
        # Start with all True
        loss_cond = df.chromosome == df.chromosome
        neut_cond = df.chromosome == df.chromosome
        gain_cond = df.chromosome == df.chromosome
        for caller in callers:
            if caller in combo:
                row[caller] = 1
                loss_cond = loss_cond & loss_conds[caller]
                neut_cond = neut_cond & neut_conds[caller]
                gain_cond = gain_cond & gain_conds[caller]
            else:
                row[caller] = 0
                loss_cond = loss_cond & ~loss_conds[caller]
                neut_cond = neut_cond & ~neut_conds[caller]
                gain_cond = gain_cond & ~gain_conds[caller]

        row['loss'] = sum(df[loss_cond].length)
        row['neutral'] = sum(df[neut_cond].length)
        row['gain'] = sum(df[gain_cond].length)

        data.append(row)

    return pd.DataFrame(data)