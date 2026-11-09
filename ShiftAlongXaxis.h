#ifndef SHIFTALONGXAXIS_H
#define SHIFTALONGXAXIS_H

#include "TGraphErrors.h" 
#include <iostream> 
using namespace std;

inline void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)
{
    // Shift original TGraphErrors along x-axis by amount determined by ’shift’.
    if(!ge){
        printf("\n WARNING: ge is NULL in ShiftAlongXaxis() !!!! \n\n");
        return;
    }
    Int_t nPoints = ge->GetN();
    Double_t x = 0.;
    Double_t y = 0.;
    for(Int_t p=0;p<nPoints;p++){
        ge->GetPoint(p,x,y);   
        x+=shift;
        ge->SetPoint(p,x,y);
    } // end of for(Int_t p=0;p<nPoints;p++)
} // end of void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)

#endif 
