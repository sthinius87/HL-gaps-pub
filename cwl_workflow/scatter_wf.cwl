#!/home/sat/anaconda3/envs/py310electronic_gaps/bin/toil-cwl-runner
cwlVersion: v1.2
class: Workflow

requirements:
  ScatterFeatureRequirement: {}

inputs:
  - id: inp_id_array
    type: int[]
    label: inp_id_array
  - id: inp_dbpath
    type: string
    label: inp_dbpath
  - id: inp_dbbasename
    type: string
    label: inp_dbbasename
  - id: inp_dbid
    type: int
    label: inp_dbid
  - id: inp_nconfs
    type: int?
    label: inp_nconfs
  - id: inp_accuracy
    type: float?
    label: inp_accuracy
  - id: inp_eltemp
    type: float?
    label: inp_eltemp
  - id: inp_method
    type: string?
    label: inp_method
  - id: inp_pathV
    type: string
  - id: inp_xtbpathV
    type: string
  - id: inp_condaprefV
    type: string
outputs:
  - id: out_results
    type: File[]
    outputSource:
      - hlgap_clt/results
#    label: out_results
steps:
  - id: hlgap_clt
    scatter: dbid
    in:
      - id: dbpath
        source: inp_dbpath
      - id: dbbasename
        source: inp_dbbasename
#      - id: dbid
#        source: inp_dbid
      - id: dbid
        source: inp_id_array
      - id: accuracy
        source: inp_accuracy
      - id: eltemp
        source: inp_eltemp
      - id: method
        source: inp_method
      - id: pathV
        source: inp_pathV
      - id: xtbpathV
        source: inp_xtbpathV
      - id: condaprefV
        source: inp_condaprefV
    out:
      - id: results
    run: ./gap.cwl
    label: step_hlgap
#    scatter: idlist
#    in:
#      idlist: inp_id_array
