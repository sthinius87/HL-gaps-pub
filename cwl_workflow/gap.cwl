class: CommandLineTool
cwlVersion: v1.2
id: hlgap
baseCommand:
  - HLgap
inputs:
  - id: dbpath
    type: string
    inputBinding:
      position: 0
      prefix: '--dbpath'
  - id: dbbasename
    type: string
    inputBinding:
      position: 0
      prefix: '--dbbasename'
  - id: dbid
    type: int
    inputBinding:
      position: 0
      prefix: '--dbid'
  - id: nconfs
    type: int?
    inputBinding:
      position: 0
      prefix: '--nconfs'
  - id: accuracy
    type: float?
    inputBinding:
      position: 0
      prefix: '--accuracy'
  - id: eltemp
    type: float?
    inputBinding:
      position: 0
      prefix: '--eltemp'
  - id: method
    type: string?
    inputBinding:
      position: 0
      prefix: '--method'
#
  - id: pathV
    type: string
  - id: xtbpathV
    type: string
  - id: condaprefV
    type: string
outputs:
  - id: results
    type: File
    outputBinding:
      glob: results*.raw
label: hlgap
requirements:
  - class: ResourceRequirement
    ramMin: 1000
    coresMin: 4
    coresMax: 12
    outdirMin: 1000
  - class: InlineJavascriptRequirement
  - class: EnvVarRequirement
    envDef:
        PATH: $(inputs.pathV)
        XTBPATH: $(inputs.xtbpathV)
        CONDA_PREFIX: $(inputs.condaprefV)