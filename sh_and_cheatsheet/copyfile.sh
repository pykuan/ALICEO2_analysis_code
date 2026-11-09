#!/bin/bash

# realdata:
(alien_cp alien:///alice/data/2023/LHC23zzh/544123/apass5/0520/o2_ctf_run00544123_orbit0181858432_tf0000163230_epn257/007/AO2D.root file:///scratch3/go52dab/O2tutorial/analysis_code/realdata/0/AO2D.root > download.log 2>&1) &
(alien_cp alien:///alice/data/2023/LHC23zzh/544123/apass5/0520/o2_ctf_run00544123_orbit0181815264_tf0000161881_epn208/001/AO2D.root file:///scratch3/go52dab/O2tutorial/analysis_code/realdata/1/AO2D.root) &
(alien_cp alien:///alice/data/2023/LHC23zzh/544123/apass5/0520/o2_ctf_run00544123_orbit0183189472_tf0000204825_epn265/001/AO2D.root file:///scratch3/go52dab/O2tutorial/analysis_code/realdata/2/AO2D.root) &
(alien_cp alien:///alice/data/2023/LHC23zzh/544123/apass5/0520/o2_ctf_run00544123_orbit0183124224_tf0000202786_epn086/001/AO2D.root file:///scratch3/go52dab/O2tutorial/analysis_code/realdata/3/AO2D.root) &
(alien_cp alien:///alice/data/2023/LHC23zzh/544123/apass5/0520/o2_ctf_run00544123_orbit0186939328_tf0000322008_epn244/001/AO2D.root file:///scratch3/go52dab/O2tutorial/analysis_code/realdata/4/AO2D.root) &
(alien_cp alien:///alice/data/2023/LHC23zzh/544123/apass5/0520/o2_ctf_run00544123_orbit0187127712_tf0000327895_epn250/005/AO2D.root file:///scratch3/go52dab/O2tutorial/analysis_code/realdata/5/AO2D.root) &
(alien_cp alien:///alice/data/2023/LHC23zzh/544123/apass5/0520/o2_ctf_run00544123_orbit0184506016_tf0000245967_epn112/002/AO2D.root file:///scratch3/go52dab/O2tutorial/analysis_code/realdata/6/AO2D.root) &
(alien_cp alien:///alice/data/2023/LHC23zzh/544123/apass5/0520/o2_ctf_run00544123_orbit0184408864_tf0000242931_epn236/001/AO2D.root file:///scratch3/go52dab/O2tutorial/analysis_code/realdata/7/AO2D.root) &

wait

echo "All jobs finished."

# MCdata: (I don't need this at this moment)