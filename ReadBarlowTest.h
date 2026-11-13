#ifndef READBARLOWTEST_H
#define READBARLOWTEST_H

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <array>

using namespace std;

// [#cent]

inline std::array<double, 8> ReadBarlowTest(const char* dataSysErr, const char* observable){

    // std::vector<std::string> observableList = {"v2{2}", "v2{4}", "v3{2}", "v3{4}", "v4{2}", "v4{4}", "NSC(2,3)", "NSC(2,4)", "NSC(2,3,4)"};
    std::vector<std::string> centBinsList = {"2.5", "7.5", "15", "25", "35", "45", "55", "65"};

    std::array<double, 8> tableBarlowTest{};
    ifstream inSysErr(dataSysErr);
    string line;
    double dCent, dSigBarlow, dSigSys;
    int index;
    bool bMatched = false;

    while (getline(inSysErr, line)) {

        if (line.find("#X,sigBarlow,sigSys", 0) == 0 || line.empty()) {
            continue;
        } // ignore comment

        else if (line.find("#end", 0) == 0) {
            bMatched = false;
            continue;
        }

        else if (line == "#:" + std::string(observable)) {
            bMatched = true;
            continue;
        } // Observable title

        if (bMatched) {
            stringstream ss(line);
            string sCent, sSigBarlow, sSigSys;

            getline(ss, sCent, ','); 
            dCent = std::stod(sCent);

            getline(ss, sSigBarlow, ','); 
            dSigBarlow = std::stod(sSigBarlow);

            getline(ss, sSigSys, ','); 
            dSigSys = std::stod(sSigSys);

            auto it = std::find(centBinsList.begin(), centBinsList.end(), sCent);
            index = std::distance(centBinsList.begin(), it);

            tableBarlowTest[index] = dSigSys;
        }
        else {
            continue;
        }

    }

    inSysErr.close();
    return tableBarlowTest;

}

#endif