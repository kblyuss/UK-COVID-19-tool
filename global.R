library(deSolve)
library(compiler)
library(Rcpp)
mm_code =
  "NumericVector my_mm(NumericMatrix m, NumericVector v){
   int nRow = m.rows();
   int nCol = m.cols();
   NumericVector ans(nRow);
   double v_j;
   for(int j = 0; j < nCol; j++){
     v_j = v[j];
     for(int i = 0; i < nRow; i++){
       ans[i] += m(i,j) * v_j;
     }
   }
   return(ans);
 }
 "
# Compiling
my_mm = cppFunction(code = mm_code)

data_regions<-read.csv("code/data_new.csv", header = TRUE)
BBC_matrix<-read.table("code/BBC_matrix_all.dat")
ph_BBC<-read.table("code/prop_hosp.dat")
pi_BBC<-read.table("code/prop_icu.dat")
pm_BBC<-read.table("code/mort.dat")

Comix_matrix<-read.table("code/mmatrix_comix_all.dat")
ph_Com<-read.table("code/prop_hosp_comix.dat")
pi_Com<-read.table("code/prop_icu_comix.dat")
pm_Com<-read.table("code/mort_comix.dat")