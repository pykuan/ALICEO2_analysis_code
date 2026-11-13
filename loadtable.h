#ifndef LOADTABLE_H
#define LOADTABLE_H

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include "TGraphErrors.h"

using namespace std;
inline vector<TGraphErrors*> loadtable(const char* csvfile)
{
    ifstream in(csvfile);

    if (!in.is_open()) {
        cerr << "ERROR: cannot open CSV file: "
             << csvfile << endl;
        return {};
    }

    vector<TGraphErrors*> graphs;
    string line, currentTitle;
    vector<double> vx, vy, ex, ey;

    auto push = [&]() {
        if (currentTitle.empty()) return;

        TGraphErrors* g = nullptr;
        if (vx.empty()) {
            g = new TGraphErrors();
        } else {
            g = new TGraphErrors(
                static_cast<int>(vx.size()),
                vx.data(), vy.data(), ex.data(), ey.data());
        }

        g->SetNameTitle(currentTitle.c_str(), currentTitle.c_str());
        graphs.push_back(g);

        currentTitle.clear();
        vx.clear(); vy.clear(); ex.clear(); ey.clear();
        // if (vx.empty()) return;
        // auto g = new TGraphErrors(vx.size(),
        //                           vx.data(), vy.data(),
        //                           ex.data(), ey.data());
        // g->SetNameTitle(currentTitle.c_str(), currentTitle.c_str());
        // graphs.push_back(g);
        // vx.clear(); vy.clear(); ex.clear(); ey.clear();
    };

    while (getline(in,line)) {
        if (line.rfind("#:",0)==0) continue;
        if (line.empty())        { push(); continue; }

        if (line.rfind("X",0)==0 || line.rfind("Centrality percentile",0)==0 || line.rfind("RUN NUMBER",0)==0) {
            stringstream hs(line);
            string dummy1, dummy2, dummy3; // CENTRALITY [PCT],CENTRALITY [PCT] LOW,X [PCT] HIGH,
            getline(hs,dummy1,',');
            getline(hs,dummy2,',');
            getline(hs,dummy3,',');
            getline(hs,currentTitle,',');
            continue;
        }

        stringstream ss(line);
        string tok;
        double cx,val,err;

        getline(ss,tok,',');  cx   = stod(tok); // x value (cent center)
        getline(ss,tok,','); // cent lower limit
        getline(ss,tok,','); // cent upper limit
        getline(ss,tok,',');  if(tok == "-") continue; else val = stod(tok); // y value (vn{k})
        getline(ss,tok,',');  err = stod(tok); // stat +
        getline(ss,tok,','); // stat -
        getline(ss,tok,',');  err += stod(tok); // sys +

        if (!std::isfinite(cx) ||
            !std::isfinite(val) ||
            !std::isfinite(err)) {
            continue;
        }

        vx.push_back(cx);
        vy.push_back(val);
        ex.push_back(0.0);
        ey.push_back(err);
    }
    push();
    return graphs;
}

#endif
