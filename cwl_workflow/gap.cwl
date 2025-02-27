cwlVersion: v1.2
class: CommandLineTool
id: hlgap

doc: |
  Calculates electronic gaps using the HLgap command-line tool.

baseCommand:
  - HLgap

inputs:
  dbpath:
    type: string
    label: Database path.
    doc: Path to the database file.
    inputBinding:
      position: 0
      prefix: "--dbpath"
  dbbasename:
    type: string
    label: Database basename.
    doc: Basename of the database file.
    inputBinding:
      position: 0
      prefix: "--dbbasename"
  dbid:
    type: int
    label: Database ID.
    doc: ID of the database entry to process.
    inputBinding:
      position: 0
      prefix: "--dbid"
  nconfs:
    type: int?
    label: Number of conformations.
    doc: Optional number of conformations to consider.
    inputBinding:
      position: 0
      prefix: "--nconfs"
  accuracy:
    type: float?
    label: Accuracy setting.
    doc: Optional accuracy setting for the calculation.
    inputBinding:
      position: 0
      prefix: "--accuracy"
  eltemp:
    type: float?
    label: Electronic temperature.
    doc: Optional electronic temperature setting.
    inputBinding:
      position: 0
      prefix: "--eltemp"
  method:
    type: string?
    label: Calculation method.
    doc: Optional method to use for the calculation.
    inputBinding:
      position: 0
      prefix: "--method"
  pathV:
    type: string
    label: Path to VASP executable.
    doc: Path to the VASP executable.
  xtbpathV:
    type: string
    label: Path to xTB executable.
    doc: Path to the xTB executable.
  condaprefV:
    type: string
    label: Conda environment prefix.
    doc: Prefix of the conda environment.

outputs:
  results:
    type: File
    label: Result file.
    doc: The output file containing the calculation results.
    outputBinding:
      glob: "results*.raw"

requirements:
  ResourceRequirement:
    ramMin: 1000
    coresMin: 4
    coresMax: 12
    outdirMin: 1000
  InlineJavascriptRequirement: {}
  EnvVarRequirement:
    envDef:
      PATH: $(inputs.pathV)
      XTBPATH: $(inputs.xtbpathV)
      CONDA_PREFIX: $(inputs.condaprefV)