#ifndef BARLOWTEST_H
#define BARLOWTEST_H

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

inline void BarlowTest(const char* dataTrial, const char* dataDefault, const char* dataSysErr){
    // const char* observableList[9] = {"v2{2}", "v2{4}", "v3{2}", "v3{4}", "v4{2}", "v4{4}", "NSC(2,3)", "NSC(2,4)", "NSC(2,3,4)"};
    ifstream inTrial(dataTrial);
    ifstream inDefault(dataDefault);
    string lineDefault, lineTrial, currentTitleDefault, currentTitleTrial;
    double dxDefault, dyDefault, derrDefault, dxTrial, dyTrial, derrTrial;
    bool bMatched = false;

    const double thresBarlow = 2.;
    double diffBarlow, errBarlow, sigBarlow, sigSys;

    ofstream outSysErr(dataSysErr);
    outSysErr << "#X,sigBarlow,sigSys" << endl;

    while (getline(inDefault, lineDefault)) {
        // printf(">> loop 1\n");
        if (lineDefault.rfind("#X,Y,staterror", 0) == 0 || lineDefault.empty()) {
            continue;
        } // ignore comment

        else if (lineDefault.rfind("#end", 0) == 0) {
            currentTitleDefault = "";
            continue;
        }

        else if (lineDefault.rfind("#", 0) == 0) {
            currentTitleDefault = lineDefault.substr(1);
            continue;
        } // Observable title

        stringstream ssD(lineDefault);
        string sxD, syD, serrD;

        getline(ssD, sxD, ','); 
        dxDefault = std::stod(sxD);

        getline(ssD, syD, ','); 
        dyDefault = std::stod(syD);

        getline(ssD, serrD, ','); 
        derrDefault = std::stod(serrD);

        // Get trial data
        while (getline(inTrial, lineTrial)) {
            // printf(">> loop 2\n");
            // printf(">> loop 2, title1 = %s, title2 = %s\n", currentTitleDefault.c_str(), currentTitleTrial.c_str());

            if (lineTrial.rfind("#X,Y,staterror", 0) == 0 || lineTrial.empty()) {
                // printf(">> loop 2, ignore comment\n");
                continue;
            } // ignore comment

            else if (lineTrial.rfind("#", 0) == 0 && lineTrial.rfind("#end", 0) != 0) {
                currentTitleTrial = lineTrial.substr(1);
                if (currentTitleDefault != currentTitleTrial) {
                    // printf(">> loop 2, # but title not matched\n");
                    continue;
                }
                else {
                    bMatched = true;
                    // printf(">> loop 2, title matched\n");
                    continue;
                }
            } // Observable title matching

            if (bMatched) {

                if (lineTrial.rfind("#end", 0) != 0) {
                    stringstream ssT(lineTrial);
                    string sxT, syT, serrT;

                    getline(ssT, sxT, ','); 
                    dxTrial = std::stod(sxT);

                    getline(ssT, syT, ','); 
                    dyTrial = std::stod(syT);

                    getline(ssT, serrT, ','); 
                    derrTrial = std::stod(serrT);

                    if (dxDefault == dxTrial) {
                        // printf("cenDef = %f, cenTrial = %f, cent matches\n", dxDefault, dxTrial);
                        diffBarlow = dyTrial - dyDefault;
                        errBarlow = std::pow(std::abs(derrTrial * derrTrial - derrDefault * derrDefault), 0.5);
                        sigBarlow = diffBarlow / errBarlow;

                        if (std::abs(sigBarlow) > thresBarlow) {
                            sigSys = diffBarlow / dyDefault; // relative systematic error
                        }
                        else {
                            sigSys = 0.;
                        }

                        if (dxTrial == 2.5) {
                            outSysErr << "#" << currentTitleTrial << endl;
                        }
                        outSysErr << dxTrial << "," << sigBarlow << "," << sigSys << endl;


                        // printf(">> loop 2, cent matched\n");
                        continue;
                    }
                    else {
                        // printf("cenDef = %f, cenTrial = %f, cent doesn't match\n", dxDefault, dxTrial);

                        // printf(">> loop 2, cent doesn't matched\n");
                        continue;
                    }
                }
                else {
                    bMatched = false;
                    if (dxDefault == 65.) {
                        outSysErr << "#end" << endl;
                    }

                    currentTitleTrial = "";
                    // printf(">> loop 2, matched, cent reach to #end\n");
                    inTrial.seekg(0);
                    break;
                }

            }
            else {

                // printf(">> loop 2, not yet match title\n");
                continue;
            }
        }

        // printf("end the second loop\n");
        dxDefault = 0.;
        dyDefault = 0.;
        derrDefault = 0.;
    }

    inTrial.close();
    inDefault.close();
    outSysErr.close();

}

#endif