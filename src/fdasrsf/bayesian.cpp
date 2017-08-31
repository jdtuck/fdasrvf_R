#include <RcppArmadillo.h>
#include <Rcpp.h>
// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

vec approx(int nd,vec xd, vec yd,int ni, vec xi);
vec findinv(mat warps, int times);
vec R_diff(vec x);

RcppExport SEXP dpcode(SEXP R_q1, SEXP R_q1L, SEXP R_q2L, SEXP R_times, SEXP R_cut)
{

  int times = as<int>(R_times);
  int cut = as<int>(R_cut);
  vec q1 = as<vec>(R_q1);
  vec q1L = as<vec>(R_q1L);
  vec q2L = as<vec>(R_q2L);

  int colnum = q1L.size();
  int rownum = colnum/times;
  int q2LLlen = (colnum-1)*times+1;
  IntegerVector tempspan1(q2LLlen);
  tempspan1.fill(0);
  tempspan1= seq_len(q2LLlen)-1;
  NumericVector temp_span2 = as<NumericVector>(tempspan1);
  float timesf = (float)(times);
  vec q2LL_time = temp_span2*(1/timesf);
  IntegerVector q2L_time1 = seq_len(colnum)-1;
  NumericVector q2L_time2 = as<NumericVector>(q2L_time1);
  vec q2LL = approx(colnum, q2L_time2, q2L, q2LLlen, q2LL_time);
  mat ID(rownum+1,colnum+1);
  ID.fill(0);
  mat S(rownum+1,colnum+1);
  S.fill(0);
  int jend,start,end,k;
  uword index;
  IntegerVector interx(times-1);
  NumericVector intery2(times-1);
  IntegerVector intery(times-1);
  IntegerVector span1 = seq_len(times-1);
  NumericVector span2 = as<NumericVector>(span1);
  vec q1x(times-1);
  vec q2y(times-1);

  for (int i = 1;i < rownum ;i++)
  {
    jend = ((i*cut+1)< colnum)?(i*cut+1):(colnum);

    for (int j=i;j< jend;j++)
    {
      start = max(IntegerVector::create(i-1,j-cut));
      end   = min(IntegerVector::create(j-1,cut*(i-1)));
      IntegerVector n = start-1+seq_len(end-start+1);
      k = LENGTH(n);
      interx = times*(i-1)+span1;
      vec Energy(k);

      for(int m = 0; m < k ; m++)
      {
        float jf = j;
        intery2 = ((jf-n[m])/times)*(span2)+n[m];
        for (int l = 0;l<(times-1);l++)
        {
          intery[l] = round(times*intery2[l]);
        }
        for (int l=0;l<(times-1);l++)
        {
          q1x[l] = q1L[interx[l]];
          q2y[l] = q2LL[intery[l]];
        }
        Energy[m] = S(i-1,n[m])+
          pow(q1[i-1]- sqrt((jf-n[m])/times)*q2L[n[m]],2)+
          pow(norm(q1x-sqrt((jf-n[m])/times)*q2y,2),2);
      }

      double min_val = Energy.min(index);
      int loc = n[index];
      S(i,j) = min_val;
      ID(i,j) = loc;
    }
  }

  int i = rownum;
  int j = colnum;
  start = max(IntegerVector::create(i-1,j-cut));
  end   = min(IntegerVector::create(j-1,cut*(i-1)));
  IntegerVector n = start-1+seq_len(end-start+1);
  k = LENGTH(n);
  interx = times*(i-1)+span1;
  span2 = as<NumericVector>(span1);
  vec Energy(k);

  for(int m = 0; m < k ; m++)
  {
    float jf = j;
    intery2 = ((jf-n[m])/times)*(span2)+n[m];
    for (int l= 0;l<(times-1);l++)
    {
      intery[l] = round(times*intery2[l]);
      if (intery[l] >= q2LLlen)
        intery[l] = q2LLlen-1;
    }
    for (int l=0;l<(times-1);l++)
    {
      q1x[l] = q1L[interx[l]];
      q2y[l] = q2LL[intery[l]];
    }
    Energy[m]=S(i-1,n[m])+
      pow(q1[i-1]-sqrt((jf-n[m])/times)*q2L[n[m]],2)+
      pow(norm(q1x-sqrt((jf-n[m])/times)*q2y,2),2);
  }

  double min_val = Energy.min(index);
  int loc = n[index];
  S(i,j) = min_val;
  ID(i,j) = loc;

  IntegerVector path(rownum);
  int count = ID(i,j);
  int oldcount;
  path(i-1) = count;

  while (count>1){
    i--;
    oldcount = count;
    count = ID(i,oldcount);
    path(i-1) = count;
  }

  IntegerVector J;
  J = path+1;

  return (List::create(Named("MatchIn2") = J,
                       Named("NDist") = S(rownum,colnum)/colnum,
                       Named("q2LL") = q2LL));
}

RcppExport SEXP simucode(SEXP R_iter, SEXP R_p, SEXP R_qt1_5, SEXP R_qt2_5,
                         SEXP R_L, SEXP R_tau, SEXP R_times, SEXP R_kappa,
                         SEXP R_alpha, SEXP R_beta, SEXP R_powera, SEXP R_dist,
                         SEXP R_dist_min, SEXP R_best_match, SEXP R_match,
                         SEXP R_thin, SEXP R_cut)
{

  int iter = as<int>(R_iter);
  int p = as<int>(R_p);
  vec match = as<vec>(R_match)-1;
  vec best_match = as<vec>(R_best_match)-1;
  int L = as<int>(R_L);
  float tau = as<float>(R_tau);
  int times = as<int>(R_times);
  float kappa = as<float>(R_kappa);
  float alpha = as<float>(R_alpha) ;
  float beta = as<float>(R_beta);
  float powera = as<float>(R_powera);
  float dist = as<float>(R_dist);
  float dist_min = as<float>(R_dist_min);
  int thin = as<int>(R_thin);
  vec qt1_5 = as<vec>(R_qt1_5);
  vec qt2_5 = as<vec>(R_qt2_5);

  IntegerVector Ixout(2*times),Ioriginal(L+1);
  IntegerVector numaccept = rep(0, L-1);
  float increment,n_dist,o_dist,adjustcon,ratio,prob,u,logpost;
  int increment_int,newj,tempnew,tempold,tempx;
  vec newmatchvec(3),oldmatchvec(3),idenmatchvec(3),pnewvec(2),poldvec(2);
  vec interynew(2*times),interyold(2*times),interx(2*times),xout(2*times),internew(2*times),interold(2*times);
  vec qt1_5_interx(2*times),qt2_5_internew(2*times),qt2_5_interold(2*times),diff_ynew(2*times),diff_yold(2*times);
  vec original(L+1),idy(p),scalevec(p),qt_5_fitted(p),kappa_collect(iter),log_collect(iter),dist_collect(iter);
  IntegerVector Irow = seq_len(p)-1;
  vec row = as<vec>(Irow);
  rowvec scale(L);
  mat match_collect(iter/thin,L+1);

  int q2LLlen = (p-1)*times+1;
  vec q2LL(q2LLlen);
  IntegerVector tempspan1(q2LLlen);
  tempspan1.fill(0);
  tempspan1 = seq_len(q2LLlen)-1;
  NumericVector temp_span2 = as<NumericVector>(tempspan1);
  float timesf = (float)(times);
  vec q2LL_time = temp_span2*(1/timesf);
  IntegerVector q2L_time1 = seq_len(p)-1;
  NumericVector q2L_time2 = as<NumericVector>(q2L_time1);
  q2LL = approx(p,q2L_time2,qt2_5,q2LLlen,q2LL_time);

  RNGScope scope;

  for (int j = 0; j < iter; j++)
  {
    for (int i = 1; i < L; i++)
    {
      if ((match[i+1]-match[i-1])>2)
      {
        increment =  as<float>(rnorm(1,0,tau));
        increment_int = round(increment);
        if (increment_int == 0) {increment_int = (increment>0)?(1):(-1);}
        newj = match[i] + increment_int;
        if((newj < match[i+1]) && (newj > match[i-1]))
        {
          newmatchvec(0) = match[i-1];
          newmatchvec(1) = newj;
          newmatchvec(2) = match[i+1];
          oldmatchvec(0) = match[i-1];
          oldmatchvec(1) = match[i];
          oldmatchvec(2) = match[i+1];
          idenmatchvec(0) = times*(i-1);
          idenmatchvec(1) = times*i;
          idenmatchvec(2) = times*(i+1);
          Ixout = seq_len(2*times)+times*(i-1)-1;
          xout = as<arma::vec>(Ixout);
          internew = approx(3,idenmatchvec,newmatchvec,2*times,xout);
          interold = approx(3,idenmatchvec,oldmatchvec,2*times,xout);
          interx = xout;
          interynew = internew;
          interynew.insert_rows(2*times,1);
          interynew[2*times] = match[i+1];
          interyold = interold;
          interyold.insert_rows(2*times,1);
          interyold[2*times] = match[i+1];
          diff_ynew = interynew.rows(1,2*times)-interynew.rows(0,2*times-1);
          diff_yold = interyold.rows(1,2*times)-interyold.rows(0,2*times-1);
          for (int ll=0;ll< (2*times); ll++)
          {
            tempx = round(interx[ll]);
            qt1_5_interx[ll] = qt1_5[tempx];
            internew[ll] = (internew[ll] > (p-1))?(p-1):(internew[ll]);
            tempnew = round(times*internew[ll]);
            qt2_5_internew[ll] = q2LL(tempnew)*sqrt(diff_ynew[ll]);
            interold[ll] = (interold[ll] > (p-1))?(p-1):(interold[ll]);
            tempold = round(times*interold[ll]);
            qt2_5_interold[ll] = q2LL[tempold]*sqrt(diff_yold[ll]);
          }
          n_dist = pow(norm(qt1_5_interx-qt2_5_internew,2),2)/p;
          o_dist = pow(norm(qt1_5_interx-qt2_5_interold,2),2)/p;
          pnewvec = (newmatchvec.rows(1,2)-newmatchvec.rows(0,1))/p;
          poldvec = (oldmatchvec.rows(1,2)-oldmatchvec.rows(0,1))/p;
          adjustcon = exp((powera-1)*(sum(log(pnewvec))-sum(log(poldvec))));
          ratio = adjustcon*exp(kappa*o_dist-kappa*n_dist);
          prob = (ratio < 1)?(ratio):(1);
          u =  as<float>(runif(1));
          if (u < prob)
          {
            match[i] = newj;
            numaccept[i-1] += 1;
          }
        }
      }
    }

    Ioriginal = (seq_len(L+1)-1)*times;
    original = as<vec>(Ioriginal);
    idy = round(approx(L+1,original,match,p,row));
    for (int ii = 0;ii<L;ii++){scale[ii] = sqrt((match[ii+1]-match[ii])/times);}
    for (int kk=0;kk<p;kk++)
    {
      idy[kk] = (idy[kk]<p)?(idy[kk]):(p-1);
      scalevec[kk] = scale[kk/times];
      qt_5_fitted[kk] = scalevec[kk]*qt2_5[idy[kk]];
    }
    dist =  pow(norm(qt1_5-qt_5_fitted,2),2)/p;
    dist_collect[j] = dist;
    if (dist < dist_min)
    {
      best_match = match;
      dist_min = dist;
    }
    if(j%thin==0){ match_collect.row(j/thin) =  match.t();}
    kappa = as<float>(rgamma(1,p/2+alpha, 1/(dist+beta)));
    kappa_collect[j] = kappa;
    logpost = (p/2+alpha)*log(kappa)-kappa*(beta+dist);
    log_collect[j] = logpost;
  }
  
  NumericVector pctaccept = as<NumericVector>(numaccept) / iter;
  vec Rr_best_match = best_match+1;
  mat Rr_match_collect = match_collect+1;

  return (List::create(Named("best_match")=Rr_best_match,
                       Named("match_collect")=Rr_match_collect,
                       Named("dist_collect")=dist_collect,
                       Named("kappa_collect")=kappa_collect,
                       Named("log_collect")=log_collect,
                       Named("dist_min")=dist_min,
                       Named("pct_accept")=pctaccept));
}

RcppExport SEXP itercode(SEXP R_iter, SEXP R_n, SEXP R_m, SEXP R_mu_5, SEXP R_match_matrix,
                         SEXP R_qt_matrix, SEXP R_qt_fitted_matrix, SEXP R_L, SEXP R_tau,
                         SEXP R_times, SEXP R_kappa, SEXP R_alpha, SEXP R_beta,
                         SEXP R_powera, SEXP R_best_vec, SEXP R_dist_vec,
                         SEXP R_best_match_matrix, SEXP R_mu_prior, SEXP R_var_const,
                         SEXP R_sumdist, SEXP R_thin, SEXP R_mu_q, SEXP R_mu_q_standard,
                         SEXP R_logmax, SEXP R_burnin, SEXP R_AVG)
{
  int iter = as<int>(R_iter);
  int n = as<int>(R_n);
  int m = as<int>(R_m);
  vec mu_5 = as<vec>(R_mu_5);
  mat match_matrix = as<mat>(R_match_matrix)-1;
  mat qt_matrix = as<mat>(R_qt_matrix);
  mat qt_fitted_matrix = as<mat>(R_qt_fitted_matrix);
  int L = as<int>(R_L);
  rowvec scale(L);
  float tau = as<float>(R_tau);
  int times = as<int>(R_times);
  float kappa = as<float>(R_kappa);
  float alpha = as<float>(R_alpha) ;
  float beta = as<float>(R_beta);
  float powera = as<float>(R_powera);
  vec best_vec = as<vec>(R_best_vec);
  vec dist_vec = as<vec>(R_dist_vec) ;
  mat best_match_matrix = as<mat>(R_best_match_matrix)-1;
  vec mu_prior = as<vec>(R_mu_prior);
  float var_const = as<float>(R_var_const);
  vec sumdist = as<vec>(R_sumdist);
  int thin = as<int>(R_thin);
  mat mu_q = as<mat>(R_mu_q);
  mat mu_q_standard = as<mat>(R_mu_q_standard);
  float logmax = as<float>(R_logmax);
  int burnin = as<int>(R_burnin);
  float AVG = as<float>(R_AVG);

  float increment,adjustcon,ratio,prob,u,n_dist,o_dist,dist,logpost,SigmaVar,rescale;
  int increment_int,newj,tempnew,tempold,tempx;
  vec newmatchvec(3),oldmatchvec(3),idenmatchvec(3),qt2_5(m),match(L+1);
  vec interynew(2*times),interyold(2*times),interx(2*times),scalevec(m),qt_5_fitted(m);
  IntegerVector Ixout(2*times),Ioriginal(L+1);
  vec xout(2*times),internew(2*times),interold(2*times),mu_5_interx(2*times);
  vec qt2_5_internew(2*times),qt2_5_interold(2*times),diff_ynew(2*times),diff_yold(2*times);
  vec pnewvec(2),poldvec(2),original(L+1),idy(m);
  IntegerVector Irow = seq_len(m)-1;
  vec row = as<vec>(Irow);
  vec karcher_res(m+1),tempkarcher_res(m+1),revscalevec(L*times),kappa_collect(iter);
  vec log_collect(iter),MAP(m),Mean(m),qt2_5_warp(m),mu_5_warp(m);
  mat tempmatch_matrix(L+1,n),bayes_warps(L+1,n);
  bayes_warps.fill(0);

  int q2LLlen = (m-1)*times+1;
  mat q2LL_mat(q2LLlen,n);
  IntegerVector tempspan1(q2LLlen);
  tempspan1.fill(0);

  tempspan1 = seq_len(q2LLlen)-1;
  NumericVector temp_span2 = as<NumericVector>(tempspan1);
  float timesf = (float)(times);
  vec q2LL_time = temp_span2*(1/timesf);
  IntegerVector q2L_time1 = seq_len(m)-1;
  NumericVector q2L_time2 = as<NumericVector>(q2L_time1);

  for (int LLL=0;LLL<n;LLL++ )
  {
    q2LL_mat.col(LLL) = approx(m,q2L_time2,qt_matrix.col(LLL),q2LLlen,q2LL_time);
  }


  RNGScope scope;

  for (int j=0; j<iter; j++)
  {
    for (int t=0; t<n; t++)
    {
      qt2_5 = qt_matrix.col(t);
      match = match_matrix.col(t);
      for (int i=1; i<L; i++)
      {
        if ((match[i+1]-match[i-1])>2)
        {
          increment =  as<float>(rnorm(1,0,tau));
          increment_int = round(increment);
          if (increment_int == 0) {break;} //{increment_int = (increment>0)?(1):(-1);}
          newj = match[i] + increment_int;
          if ((newj < match[i+1]) && (newj > match[i-1]))
          {
            newmatchvec(0) = match[i-1];
            newmatchvec(1) = newj;
            newmatchvec(2) = match[i+1];
            oldmatchvec(0) = match[i-1];
            oldmatchvec(1) = match[i];
            oldmatchvec(2) = match[i+1];
            idenmatchvec(0) = times*(i-1);
            idenmatchvec(1) = times*i;
            idenmatchvec(2) = times*(i+1);
            Ixout = seq_len(2*times)+times*(i-1)-1;
            xout = as<arma::vec>(Ixout);
            internew = approx(3,idenmatchvec,newmatchvec,2*times,xout);
            interold = approx(3,idenmatchvec,oldmatchvec,2*times,xout);
            interx = xout;
            interynew = internew;
            interynew.insert_rows(2*times,1);
            interynew[2*times] = match[i+1];
            interyold = interold;
            interyold.insert_rows(2*times,1);
            interyold[2*times] = match[i+1];
            diff_ynew = interynew.rows(1,2*times)-interynew.rows(0,2*times-1);
            diff_yold = interyold.rows(1,2*times)-interyold.rows(0,2*times-1);
            for (int ll=0;ll< (2*times); ll++)
            {
              tempx = round(interx[ll]);
              mu_5_interx[ll] = mu_5[tempx];
              internew[ll] = (internew[ll] > (m-1))?(m-1):(internew[ll]);
              tempnew = round(times*internew[ll]);
              tempnew = (tempnew>q2LLlen-1)?(q2LLlen-1):(tempnew);
              qt2_5_internew[ll] = q2LL_mat(tempnew,t)*sqrt(diff_ynew[ll]);
              interold[ll] = (interold[ll] > (m-1))?(m-1):(interold[ll]);
              tempold = round(times*interold[ll]);
              tempold = (tempold>q2LLlen-1)?(q2LLlen-1):(tempold);
              qt2_5_interold[ll] = q2LL_mat(tempold,t)*sqrt(diff_yold[ll]);
            }
            n_dist = pow(norm(mu_5_interx-qt2_5_internew,2),2)/m;
            o_dist = pow(norm(mu_5_interx-qt2_5_interold,2),2)/m;
            pnewvec = (newmatchvec.rows(1,2)-newmatchvec.rows(0,1))/m;
            poldvec = (oldmatchvec.rows(1,2)-oldmatchvec.rows(0,1))/m;
            adjustcon = exp((powera-1)*(sum(log(pnewvec))-sum(log(poldvec))));
            ratio = adjustcon*exp(kappa*o_dist-kappa*n_dist);
            prob = (ratio < 1)?(ratio):(1);
            u =  as<float>(runif(1));
            match[i] = (u < prob)?(newj):(match[i]);
          }
        }
      }
      Ioriginal = (seq_len(L+1)-1)*times;
      original = as<vec>(Ioriginal);
      idy = approx(L+1,original,match,m,row);
      qt2_5_warp = approx(m,q2L_time2,qt2_5,m,idy);

      for (int ii = 0;ii<L;ii++){scale[ii] = sqrt((match[ii+1]-match[ii])/times);}
      for (int kk=0;kk<m;kk++)
      {
        scalevec[kk] = scale[kk/times];
        qt_5_fitted[kk] = scalevec[kk]*qt2_5_warp[kk];
      }
      dist =  pow(norm(mu_5-qt_5_fitted,2),2)/m;
      match_matrix.col(t) = match;
      dist_vec[t] = dist;
      qt_fitted_matrix.col(t) = qt_5_fitted;
      if ((j/thin)>=burnin) {bayes_warps.col(t) = bayes_warps.col(t) + match/AVG;}
    }
    sumdist[j] = sum(dist_vec);
    if (sumdist[j] < sum(best_vec))
    {
      best_vec = dist_vec;
      best_match_matrix = match_matrix;
    }
    if(j%thin==0)
    {
      mu_q.row(j/thin) = conv_to<rowvec>::from(mu_5);
      karcher_res = findinv(match_matrix,times);
      for (int jj=0;jj<m;jj++){revscalevec[jj] = sqrt(karcher_res(jj+1)-karcher_res(jj));}
      mu_5_warp =  approx(m,q2L_time2,mu_5,m,karcher_res);
      for (int qq=0;qq<m;qq++){mu_q_standard(j/thin,qq) = revscalevec[qq]*mu_5_warp[qq];}
    }
    kappa = as<float>(rgamma(1,n*m/2+alpha, 1/(sum(dist_vec)+beta)));
    kappa_collect[j] = kappa;
    logpost = (n*m/2+alpha)*log(kappa)-kappa*(beta+sum(dist_vec));
    log_collect[j] = logpost;
    if(logpost > logmax){logmax = logpost;MAP = mu_5;}
    Mean = (var_const*m)/(m+2*kappa*n*var_const)*(2*kappa*n*mean(qt_fitted_matrix,1)/m +mu_prior/var_const);
    SigmaVar = (var_const*m)/(m+2*kappa*n*var_const);
    mu_5 = Mean + as<vec>(rnorm(m,0,sqrt(SigmaVar)));
    rescale = sqrt(m/sum(square(mu_5)));
    mu_5 = rescale*mu_5;
  }

  mat Rr_best_match_matrix = best_match_matrix+1;
  mat R_bayes_warps = bayes_warps+1;

  return List::create(Named("mu.q.standard") = mu_q_standard,
                      Named("mu.q") = mu_q,
                      Named("best_match.matrix") = Rr_best_match_matrix,
                      Named("kappafamily") = kappa_collect,
                      Named("sumdist") = sumdist,
                      Named("log.posterior") = log_collect,
                      Named("dist") = sum(best_vec), Named("MAP") = MAP,
                      Named("bayes_warps") = R_bayes_warps);
}

vec approx(int nd,vec xd, vec yd, int ni, vec xi)
{
  int i,k;
  double t;
  vec yi(ni);
  for ( i = 0; i < ni; i++ )
  {
    if ( xi(i) <= xd(0) )
    {
      t = ( xi(i) - xd(0) ) / ( xd(1) - xd(0) );
      yi(i) = ( 1.0 - t ) * yd(0) + t * yd(1);
    }
    else if ( xd(nd-1) <= xi(i) )
    {
      t = ( xi(i) - xd(nd-2) ) / ( xd(nd-1) - xd(nd-2) );
      yi(i) = ( 1.0 - t ) * yd(nd-2) + t * yd(nd-1);
    }
    else
    {
      for ( k = 1; k < nd; k++ )
      {
        if ( xd(k-1) <= xi(i) && xi(i) <= xd(k) )
        {
          t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) );
          yi(i) = ( 1.0 - t ) * yd(k-1) + t * yd(k);
          break;
        }
      }
    }
  }

  return yi;

}

vec findinv(mat warps, int times)
{

  int m = warps.n_rows;
  int n = warps.n_cols;
  mat psi_m(m-1,n);
  psi_m.zeros();

  for(int j = 0; j<n; j++)
  {
    for (int hh = 0; hh< m-1 ; hh++)
    {
      psi_m(hh,j) = sqrt( (warps(hh+1,j)-warps(hh,j))/times);
    }
  }

  vec w = mean(psi_m,1);
  vec mupsi = w/as_scalar(sqrt(sum(pow(w,2)/(m-1))));
  mat v_m(m-1,n);
  v_m.zeros();
  vec mupsi_update(m-1);
  float check = 1;

  while (check > 0.01)
  {
    for (int i = 0;i<n;i++)
    {
      double theta = acos(sum(mupsi.t()*psi_m.col(i)/(m-1)));
      v_m.col(i) = theta/sin(theta)*(psi_m.col(i)-cos(theta)*mupsi);
    }
    vec vbar = mean(v_m,1);
    check = norm(vbar,2)/sqrt(m-1);

    if (check > 0)
    {
      mupsi_update = cos(0.01*norm(vbar,2)/sqrt(m-1))*mupsi
      +sin(0.01*norm(vbar,2)/sqrt(m-1))*vbar/(norm(vbar,2)/sqrt(m-1));
    }
    else
    {
      mupsi_update = cos(0.01*norm(vbar,2)/sqrt(m-1))*mupsi;
    }
  }
  vec karcher_s = cumsum(pow(mupsi_update,2))*times;
  karcher_s.insert_rows(0,1);
  IntegerVector Irow = (seq_len(m)-1)*times;
  vec row = as<vec>(Irow);
  IntegerVector Iout = seq_len((m-1)*times+1)-1;
  vec out = as<vec>(Iout);
  vec invidy = approx(m,karcher_s,row,(m-1)*times+1,out);

  return invidy;

}

vec R_diff(vec x)
{
  int n = x.n_elem;
  vec x_diff(n-1);

  for (int i=0;i<(n-1);i++)
  {
    x_diff[i] = x[i+1]-x[i];
  }

  return x_diff;

}
