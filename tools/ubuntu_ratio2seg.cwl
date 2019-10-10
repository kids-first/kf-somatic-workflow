cwlVersion: v1.0
class: CommandLineTool
id: ubuntu_ratio2seg
requirements:
  - class: DockerRequirement
    dockerPull: 'kfdrc/python:2.7.13'
  - class: ResourceRequirement
    ramMin: 1000
    coresMin: 1
  - class: InlineJavascriptRequirement
baseCommand: [python, -c]
arguments:
  - position: 0
    valueFrom: |

        import math
        
        fai = open("$(inputs.reference_fai.path)")
        h = {}
        for line in fai:
            f = line.split("\t")
            if f == "chrM":
                break
            f[0] = f[0].replace("chr","")
            h[f[0]] = f[1]
        fai.close()

        smp = "$(inputs.sample_name)"

        ratio_file = open("$(inputs.ctrlfreec_ratio.path)")
        out = open("$(inputs.output_basename).controlfreec.seg", "w")
        out.write("ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n") 
        head = next(ratio_file)
        count = 0
        for line in ratio_file:
            data = line.rstrip("\n").split("\t")
            (chrom, pos, ratio, meanRatio) = (data[0], data[1], data[2], data[3])
            if meanRatio == "-1":
                continue
            count += 1
            if count == 1:
                start = pos
                seg_ratio = meanRatio
                on_chr = chrom
            else:
                if chrom != on_chr:
                    out.write("\t".join((smp, "chr" + on_chr, start, h[on_chr], str(count))) + "\t")
                    if seg_ratio != "0":
                        out.write(str(math.log(float(seg_ratio), 2)) + "\n")
                    else:
                        out.write(str(math.log(float(seg_ratio) + 1, 2)) + "\n")
                    start = pos
                    seg_ratio = meanRatio
                    on_chr = chrom
                    count = 1
                elif meanRatio != seg_ratio:
                    out.write("\t".join((smp, "chr" + chrom, start, str(int(pos)-1), str(count))) + "\t")
                    if seg_ratio != "0":
                        out.write(str(math.log(float(seg_ratio), 2)) + "\n")
                    else:
                        out.write(str(math.log(float(seg_ratio) + 1, 2)) + "\n")
                    start = pos
                    seg_ratio = meanRatio
                    count = 1
        ratio_file.close()
        out.close()
inputs:
  reference_fai: File
  ctrlfreec_ratio: File
  sample_name: string
  output_basename: string

outputs:
  ctrlfreec_ratio2seg:
    type: File
    outputBinding:
      glob: '*.seg'
