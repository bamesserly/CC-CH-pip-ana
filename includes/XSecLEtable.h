//For the tables, the first raw is the xsec CV (/x10-42cm^2/MeV/nucleon), 
//the second raw is the absolute total error.
//The "Old" version is obtained from the paper https://arxiv.org/pdf/1406.6415.pdf
//The "New" version is obtained from the updated table 
//https://minerva.fnal.gov/wp-content/uploads/2017/03/Updated_1pi_data.pdf

#ifndef XSecLEtable_h
#define XSecLEtable_h

double tpiOld[2][7] = {{11.3, 11.6, 10.7, 8.50, 7.60, 6.60, 3.80},
		       {3.00, 2.50, 2.00, 1.50, 1.40, 1.10, 0.80}};
double tpiNew[2][7] = {{12.0, 13.0, 12.0, 9.00, 8.00, 7.00, 4.00},
		       {3.12, 2.73, 2.16, 1.53, 1.44, 1.19, 0.88}};
double thetapiOld[2][13] = {{18.3, 28.7, 30.5, 38.7, 35.4, 29.1, 21.3, 19.8, 15.5, 9., 7.1, 5.4, 3.3},
                            {  4.,  5.7,   6.,  7.7,  7.4,  6.1,  4.5,  4.0,  2.9,1.9, 1.4, 1.1, 0.7}};
double thetapiNew[2][13] = {{ 12.,  23.,  28.,  38.,  36.,  30.,  23.,  22.,  17.,10., 8. ,  6.,  4.},
			    {2.76, 4.83,  5.6,  7.6,  7.2,   6., 4.83,  4.4, 3.23,2.1,1.52,1.14,0.84}};

double GetXSecLECV(std::string var, int bin){
  if {var == "tpi_old"} return tpiOld[0][bin];
  else if (var == "tpi_new") return tpiNew[0][bin];
  else if (var == "thetapi_old") return thetapiOld[0][bin];
  else if (var == "thetapi_new") return thetapiNew[0][bin];
}

double GetXSecLEerror(std::string var, int bin){
  if {var == "tpi_old"} return tpiOld[1][bin];
  else if (var == "tpi_new") return tpiNew[1][bin];
  else if (var == "thetapi_old") return thetapiOld[1][bin];
  else if (var == "thetapi_new") return thetapiNew[1][bin];
}




#endif 
