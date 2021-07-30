/* Copyright (c) David Dalrymple 2011 */
#define TAU 6.283185307179586476925287 // also known as "two pi" to the unenlightened

double finv(double t) {
  return (t>(6.0/29.0))?(t*t*t):(3*(6.0/29.0)*(6.0/29.0)*(t-4.0/29.0));
}

/* Convert from L*a*b* doubles to XYZ doubles
 * Formulas drawn from http://en.wikipedia.org/wiki/Lab_color_space
 */
void lab2xyz(double* x, double* y, double* z, double l, double a, double b) {
  double sl = (l+0.16)/1.16;
  double ill[3] = {0.9643,1.00,0.8251}; //D50
  *y = ill[1] * finv(sl);
  *x = ill[0] * finv(sl + (a/5.0));
  *z = ill[2] * finv(sl - (b/2.0));
}

double correct(double cl) {
  double a = 0.055;
  return (cl<=0.0031308)?(12.92*cl):((1+a)*pow(cl,1/2.4)-a);
}

/* Convert from XYZ doubles to sRGB bytes
 * Formulas drawn from http://en.wikipedia.org/wiki/Srgb
 * Returns whether the value was out of the RGB gamut
 */
bool xyz2rgb(unsigned char* r, unsigned char* g, unsigned char* b, double x, double y, double z, bool clipOutOfGamut=false) {
  double rl =  3.2406*x - 1.5372*y - 0.4986*z;
  double gl = -0.9689*x + 1.8758*y + 0.0415*z;
  double bl =  0.0557*x - 0.2040*y + 1.0570*z;
  bool needsClip = (rl < 0.0 || rl > 1.0 || gl < 0.0 || gl > 1.0 || bl < 0.0 || bl > 1.0);
  if(needsClip) {
    rl = (rl<0.0)?0.0:((rl>1.0)?1.0:rl);
    gl = (gl<0.0)?0.0:((gl>1.0)?1.0:gl);
    bl = (bl<0.0)?0.0:((bl>1.0)?1.0:bl);
  }
  //Uncomment the below to detect clipping by making clipped zones white
  if(clipOutOfGamut && needsClip) {rl=1.0;gl=1.0; bl=1.0;}
  *r = (unsigned char)(255.0*correct(rl));
  *g = (unsigned char)(255.0*correct(gl));
  *b = (unsigned char)(255.0*correct(bl));

  return needsClip;
}

/* Convert from LAB doubles to sRGB bytes 
 * (just composing the above transforms)
 * Returns whether the value was out of the RGB gamut
 */
bool lab2rgb(unsigned char* R, unsigned char* G, unsigned char* B, double l, double a, double b, bool clipOutOfGamut=false) {
  double x,y,z;
  lab2xyz(&x,&y,&z,l,a,b);
  return xyz2rgb(R,G,B,x,y,z, clipOutOfGamut);
}

/* Convert from a qualitative parameter c and a quantitative parameter l to a 24-bit pixel
 * These formulas were invented by me to obtain maximum contrast without going out of gamut
 * if the parameters are in the range 0-1
 */
void cl2pix(void* rgb, double c, double l) {
  unsigned char* ptr = (unsigned char*)rgb;
  double L = l*0.61+0.09; //L of L*a*b*
  double angle = TAU/6.0-c*TAU;
  double r = l*0.311+0.125; //~chroma
  double a = sin(angle)*r;
  double b = cos(angle)*r;
  lab2rgb(ptr,ptr+1,ptr+2,L,a,b);
}
