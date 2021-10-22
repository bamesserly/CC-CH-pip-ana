#ifndef MYPLOTSTYLE_H
#define MYPLOTSTYLE_H

#include "TStyle.h"
#include "TColor.h"
#include "lab2rgb.h"

#include "TList.h"
#include "TPad.h"
#include "TH1.h"
#include "TCanvas.h"

#include <iostream>

//======================================================================
void setLabPalette(double lStart=0.1, double lEnd=0.9)
{
  // A colour palette based on the La*b* colour space, and stolen from http://davidad.net/colorviz/ 
  const int nCols=20;
  gStyle->SetNumberContours(999);
  Double_t stops[nCols];
  Double_t reds[nCols];
  Double_t greens[nCols];
  Double_t blues[nCols];

  for(int i=0; i<nCols; ++i){
    float frac=float(i)/nCols;
    stops[i]=frac;

    const double l = lStart+frac*(lEnd-lStart);
    double r=l*(1-0.1)+0.1;
    const double angle = frac*TAU + TAU/2;
    bool wasClipped=true;

    unsigned char red, green, blue;
    while(wasClipped && r>0){
      const double a = sin(angle)*r;
      const double b = cos(angle)*r;
      wasClipped=lab2rgb(&red, &green, &blue,
                         l, a, b, 
                         true);
      r-=0.01;
    }
    reds[i]=red/255.;
    greens[i]=green/255.;
    blues[i]=blue/255.;
  }

  TColor::CreateGradientColorTable(nCols, stops, reds, greens, blues, 999);
}

//======================================================================
void setRedPalette()
{
  const int NRGBs = 9;
  const int NCont = 99;
  static bool initialized=false;
  static int colors[99];

  if(!initialized){
    // White -> red
    Double_t stops[NRGBs] = { 0.00, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000};
    Double_t red[NRGBs]   = { 1.00, 1.00, 0.99, 0.99, 0.98, 0.94, 0.80, 0.65, 0.40 };
    Double_t green[NRGBs] = { 0.96, 0.88, 0.73, 0.57, 0.42, 0.23, 0.09, 0.06, 0.00 };
    Double_t blue[NRGBs]  = { 0.94, 0.82, 0.63, 0.45, 0.29, 0.17, 0.11, 0.08, 0.05 };
    int colmin=TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    for(int i=0; i<NCont; ++i) colors[i]=colmin+i;

    initialized=true;
  }
  gStyle->SetNumberContours(NCont);
  gStyle->SetPalette(NCont, colors);

}

//======================================================================
void setCorrelationPalette(double whiteFrac=0.5)
{
  // A colour palette that goes blue->white->red, useful for
  // correlation matrices
  const int NRGBs = 3;
  const int NCont = 99;
  static bool initialized=false;
  static int colors[99];

  if(!initialized){
    gStyle->SetNumberContours(NCont);
    Double_t stops[NRGBs] = { 0.00, whiteFrac, 1.00};
    Double_t red[NRGBs]   = { 0.00, 1.00,      1.00};
    Double_t green[NRGBs] = { 0.00, 1.00,      0.00};
    Double_t blue[NRGBs]  = { 1.00, 1.00,      0.00};
    int colmin=TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    std::cout << "colmin is " << colmin << std::endl;
    for(int i=0; i<NCont; ++i) colors[i]=colmin+i;

    initialized=true;
  }
  gStyle->SetNumberContours(NCont);
  gStyle->SetPalette(NCont, colors);
}


//======================================================================
void setCoolWarmPalette()
{
  // The "coolwarm" palette from matplotlib
  const int NRGBs = 21;
  const int NCont = 99;
  static bool initialized=false;
  static int colors[99];

  if(!initialized){
    gStyle->SetNumberContours(NCont);
    Double_t stops[NRGBs] = { 0.  ,  0.05,  0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,
                              0.45,  0.5 ,  0.55,  0.6 ,  0.65,  0.7 ,  0.75,  0.8 ,  0.85,
                              0.9 ,  0.95,  1.  };
    Double_t red[NRGBs]   = {0.23 ,  0.285,  0.348,  0.415,  0.484,  0.554,  0.619,  0.688,
                             0.754,  0.814,  0.867,  0.913,  0.947,  0.966,  0.968,  0.957,
                             0.932,  0.892,  0.839,  0.774,  0.706};
    Double_t green[NRGBs] = {0.299,  0.38 ,  0.466,  0.547,  0.622,  0.69 ,  0.744,  0.793,
                             0.83 ,  0.854,  0.864,  0.837,  0.795,  0.74 ,  0.674,  0.598,
                             0.519,  0.425,  0.322,  0.2  ,  0.016};
    Double_t blue[NRGBs]  = {0.754,  0.823,  0.888,  0.939,  0.975,  0.996,  0.999,  0.988,
                             0.961,  0.918,  0.863,  0.795,  0.717,  0.637,  0.557,  0.477,
                             0.406,  0.333,  0.265,  0.203,  0.15};
    int colmin=TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    std::cout << "colmin is " << colmin << std::endl;
    for(int i=0; i<NCont; ++i) colors[i]=colmin+i;

    initialized=true;
  }
  gStyle->SetNumberContours(NCont);
  gStyle->SetPalette(NCont, colors);

}

//======================================================================
void setRainbowToWhitePalette()
{
  // Matt Strait's colour palette that fades out to white at zero
  const int NRGBs = 7, NCont = 999;
  gStyle->SetNumberContours(NCont);
  Double_t stops[NRGBs] = { 0.00, 0.05, 0.23, 0.45, 0.60, 0.85, 1.00 };
  Double_t red[NRGBs]   = { 1.00, 0.00, 0.00, 0.00, 1.00, 1.00, 0.33 };
  Double_t green[NRGBs] = { 1.00, 1.00, 0.30, 0.40, 1.00, 0.00, 0.00 };
  Double_t blue[NRGBs]  = { 1.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
}

//======================================================================
void setBlackbodyPalette()
{
  // Available as gStyle->SetPalette(56) in sufficiently recent versions of ROOT, but not ours
  const int nRGBs = 5;
  const int NCont = 99;
  static bool initialized=false;
  static int colors[99];

  if(!initialized){
    gStyle->SetNumberContours(NCont);
    Double_t stops[nRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00};
    Double_t red[nRGBs]   = { 1.00, 1.00, 1.00, 0.50, 0.00};
    Double_t green[nRGBs] = { 1.00, 1.00, 0.55, 0.00, 0.00};
    Double_t blue[nRGBs]  = { 1.00, 0.00, 0.00, 0.00, 0.00};
    int colmin=TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, NCont);
    for(int i=0; i<NCont; ++i) colors[i]=colmin+i;

    initialized=true;
  }
  gStyle->SetNumberContours(NCont);
  gStyle->SetPalette(NCont, colors);
}

//======================================================================
void setBlackbodyPaletteInverted()
{
  // Available as gStyle->SetPalette(56) in sufficiently recent versions of ROOT, but not ours
  const int nRGBs = 5;
  const int NCont = 99;
  static bool initialized=false;
  static int colors[99];

  if(!initialized){
    gStyle->SetNumberContours(NCont);
    Double_t stops[nRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00};
    Double_t red[nRGBs]   = { 0.00, 0.50, 1.00, 1.00, 1.00};
    Double_t green[nRGBs] = { 0.00, 0.00, 0.55, 1.00, 1.00};
    Double_t blue[nRGBs]  = { 0.00, 0.00, 0.00, 0.00, 1.00};
    int colmin=TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, NCont);
    for(int i=0; i<NCont; ++i) colors[i]=colmin+i;

    initialized=true;
  }
  gStyle->SetNumberContours(NCont);
  gStyle->SetPalette(NCont, colors);
}

//======================================================================
void myPlotStyle()
{
  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadRightMargin(0.15);
  //gStyle->SetPadTopMargin(0.05);
  //gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadLeftMargin(0.15);

  // use large fonts

  gStyle->SetLabelFont(42);
  gStyle->SetTitleFont(42);

  gStyle->SetTextSize(0.08);

  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.06,"z");
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleX(0.25);
  gStyle->SetTitleY(0.95);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleOffset(0.9, "Y");

  // use bold lines and markers
  gStyle->SetMarkerStyle(1);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  //gStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(1);
  //gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  gStyle->SetTickLength(0.01, "Y");
  gStyle->SetTickLength(0.02, "X");

  gStyle->SetNdivisions(505, "XYZ");
  gStyle->SetStripDecimals(false);
}
#endif
