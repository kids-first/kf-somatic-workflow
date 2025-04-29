cwlVersion: v1.2
class: CommandLineTool
id: clt_flatten_file_list
requirements:
  - class: InlineJavascriptRequirement
baseCommand: [echo, done]
inputs:
  input_list: { type: 'File[]?' }
outputs:
  output:
    type: File[]?
    outputBinding:
      outputEval: |
        ${
          var flatten = function flatten(ary) {
              var ret = [];
              if (ary != null){
                for(var i = 0; i < ary.length; i++) {
                    if(Array.isArray(ary[i])) {
                        ret = ret.concat(flatten(ary[i]));
                    } else {
                        ret.push(ary[i]);
                    }
                }
              }
              return ret;
          }
          var flatin = flatten(inputs.input_list);
          return flatin;
        }
