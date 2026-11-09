#include "BarlowTest.h"
using namespace std;
int BarlowTestC() {
    BarlowTest("result/barlow/result_phipt_realdata.csv", "result/barlow/result_all_realdata.csv", "result/barlow/barlow_phiptVSall_real.csv");
    BarlowTest("result/barlow/result_phi_realdata.csv", "result/barlow/result_all_realdata.csv", "result/barlow/barlow_phiVSall_realdata.csv");
    BarlowTest("result/barlow/result_phipt_MCdata.csv", "result/barlow/result_all_MCdata.csv", "result/barlow/barlow_phiptVSall_MCdata.csv");
    BarlowTest("result/barlow/result_pt_MCdata.csv", "result/barlow/result_all_MCdata.csv", "result/barlow/barlow_ptVSall_MCdata.csv");
    return 0;
}