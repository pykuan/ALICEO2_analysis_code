#!/bin/bash
cp /scratch3/go52dab/O2tutorial/analysis_code/multiparticleCumulants.cxx /scratch3/go52dab/aliceO2-2/O2Physics/PWGCF/MultiparticleCorrelations/Tasks/
cp /scratch3/go52dab/aliceO2-2/O2Physics/PWGCF/MultiparticleCorrelations/Tasks/multiparticleCumulants.cxx /scratch3/go52dab/aliceO2-fork/O2Physics/PWGCF/MultiparticleCorrelations/Tasks/multiparticleCumulants.cxx
cp /scratch3/go52dab/O2tutorial/analysis_code/multiparticleCumulants.cxx /scratch3/go52dab/aliceO2-fork/O2Physics/PWGCF/MultiparticleCorrelations/Tasks/multiparticleCumulants.cxx

cd /scratch3/go52dab/aliceO2-2
aliBuild build O2Physics -j 8 --defaults o2 --debug

#cd /scratch3/go52dab/aliceO2-2/O2Physics
#clang-format -i /scratch3/go52dab/O2tutorial/analysis_code/multiparticleCumulants.cxx
#python3 /scratch3/go52dab/aliceO2-2/O2Physics/Scripts/o2_linter.py /scratch3/go52dab/O2tutorial/analysis_code/multiparticleCumulants.cxx | tee /scratch3/go52dab/O2tutorial/analysis_code/log_linter.log

# git remote add pykuan git@github.com:pykuan/O2Physics
# git remote -v
# alternative:
# sync fork repo and discard previous commit
# git reset --hard upstream/master



# cp /scratch3/go52dab/O2tutorial/analysis_code/multiparticleCumulants.cxx /scratch3/go52dab/aliceO2-fork/O2Physics/PWGCF/MultiparticleCorrelations/Tasks/multiparticleCumulants.cxx
# git status
# git add .
# git commit -m "..."
# git push pykuan

