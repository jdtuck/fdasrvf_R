#include "bayesian.h"

Rcpp::List dpcode(arma::vec q1, arma::vec q1L, arma::vec q2L, int times, int cut) {
  int colnum = q1L.size();
  int rownum = colnum/times;
  int q2LLlen = (colnum-1)*times+1;
  Rcpp::IntegerVector tempspan1(q2LLlen);
  tempspan1.fill(0);
  tempspan1 = Rcpp::seq_len(q2LLlen)-1;
  Rcpp::NumericVector temp_span2 = Rcpp::as<Rcpp::NumericVector>(tempspan1);
  float timesf = (float)(times);
  arma::vec q2LL_time = temp_span2*(1/timesf);
  Rcpp::IntegerVector q2L_time1 = Rcpp::seq_len(colnum)-1;
  Rcpp::NumericVector q2L_time2 = Rcpp::as<Rcpp::NumericVector>(q2L_time1);
  arma::vec q2LL = approx(colnum, q2L_time2, q2L, q2LLlen, q2LL_time);
  arma::mat ID(rownum+1,colnum+1);
  ID.fill(0);
  arma::mat S(rownum+1,colnum+1);
  S.fill(0);
  int jend,start,end,k;
  arma::uword index;
  Rcpp::IntegerVector interx(times-1);
  Rcpp::NumericVector intery2(times-1);
  Rcpp::IntegerVector intery(times-1);
  Rcpp::IntegerVector span1 = Rcpp::seq_len(times-1);
  Rcpp::NumericVector span2 = Rcpp::as<Rcpp::NumericVector>(span1);
  arma::vec q1x(times-1);
  arma::vec q2y(times-1);

  for (int i = 1;i < rownum ;i++)
  {
    jend = ((i*cut+1)< colnum)?(i*cut+1):(colnum);

    for (int j=i;j< jend;j++)
    {
      start = max(Rcpp::IntegerVector::create(i-1,j-cut));
      end   = min(Rcpp::IntegerVector::create(j-1,cut*(i-1)));
      Rcpp::IntegerVector n = start - 1 + Rcpp::seq_len(end-start+1);
      k = LENGTH(n);
      interx = times*(i-1)+span1;
      arma::vec Energy(k);

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

      double min_val = Energy.min();
      index = Energy.index_min();
      int loc = n[index];
      S(i,j) = min_val;
      ID(i,j) = loc;
    }
  }

  int i = rownum;
  int j = colnum;
  start = max(Rcpp::IntegerVector::create(i-1,j-cut));
  end   = min(Rcpp::IntegerVector::create(j-1,cut*(i-1)));
  Rcpp::IntegerVector n = start - 1 + Rcpp::seq_len(end-start+1);
  k = LENGTH(n);
  interx = times*(i-1)+span1;
  span2 = Rcpp::as<Rcpp::NumericVector>(span1);
  arma::vec Energy(k);

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

  double min_val = Energy.min();
  index = Energy.index_min();
  int loc = n[index];
  S(i,j) = min_val;
  ID(i,j) = loc;

  Rcpp::IntegerVector path(rownum);
  int count = ID(i,j);
  int oldcount;
  path(i-1) = count;

  while (count>1){
    i--;
    oldcount = count;
    count = ID(i,oldcount);
    path(i-1) = count;
  }

  Rcpp::IntegerVector J;
  J = path+1;

  return (Rcpp::List::create(Rcpp::Named("MatchIn2") = J,
                             Rcpp::Named("NDist") = S(rownum,colnum)/colnum,
                             Rcpp::Named("q2LL") = q2LL));
}

Rcpp::List simucode(int iter, int p, arma::vec qt1_5, arma::vec qt2_5, int L,
                    float tau, int times, float kappa, float alpha, float beta,
                    float powera, float dist, float dist_min,
                    arma::vec best_match, arma::vec match, int thin, int cut)
{
  match = match - 1;
  best_match = best_match - 1;

  Rcpp::IntegerVector Ixout(2*times),Ioriginal(L+1);
  Rcpp::IntegerVector numaccept = Rcpp::rep(0, L-1);
  float increment,n_dist,o_dist,adjustcon,ratio,prob,u,logpost;
  int increment_int,newj,tempnew,tempold,tempx;
  arma::vec newmatchvec(3),oldmatchvec(3),idenmatchvec(3),pnewvec(2),poldvec(2);
  arma::vec interynew(2*times),interyold(2*times),interx(2*times),xout(2*times),internew(2*times),interold(2*times);
  arma::vec qt1_5_interx(2*times),qt2_5_internew(2*times),qt2_5_interold(2*times),diff_ynew(2*times),diff_yold(2*times);
  arma::vec original(L+1),idy(p),scalevec(p),qt_5_fitted(p),kappa_collect(iter),log_collect(iter),dist_collect(iter);
  Rcpp::IntegerVector Irow = Rcpp::seq_len(p)-1;
  arma::vec row = Rcpp::as<arma::vec>(Irow);
  arma::rowvec scale(L);
  arma::mat match_collect(iter/thin,L+1);

  int q2LLlen = (p-1)*times+1;
  arma::vec q2LL(q2LLlen);
  Rcpp::IntegerVector tempspan1(q2LLlen);
  tempspan1.fill(0);
  tempspan1 = Rcpp::seq_len(q2LLlen)-1;
  Rcpp::NumericVector temp_span2 = Rcpp::as<Rcpp::NumericVector>(tempspan1);
  float timesf = (float)(times);
  arma::vec q2LL_time = temp_span2*(1/timesf);
  Rcpp::IntegerVector q2L_time1 = Rcpp::seq_len(p)-1;
  Rcpp::NumericVector q2L_time2 = Rcpp::as<Rcpp::NumericVector>(q2L_time1);
  q2LL = approx(p,q2L_time2,qt2_5,q2LLlen,q2LL_time);

  Rcpp::RNGScope scope;

  for (int j = 0; j < iter; j++)
  {
    for (int i = 1; i < L; i++)
    {
      if ((match[i+1]-match[i-1])>2)
      {
        increment =  Rcpp::as<float>(Rcpp::rnorm(1,0,tau));
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
          Ixout = Rcpp::seq_len(2*times)+times*(i-1)-1;
          xout = Rcpp::as<arma::vec>(Ixout);
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
          u =  Rcpp::as<float>(Rcpp::runif(1));
          if (u < prob)
          {
            match[i] = newj;
            numaccept[i-1] += 1;
          }
        }
      }
    }

    Ioriginal = (Rcpp::seq_len(L+1)-1)*times;
    original = Rcpp::as<arma::vec>(Ioriginal);
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
    kappa = Rcpp::as<float>(Rcpp::rgamma(1,p/2+alpha, 1/(dist+beta)));
    kappa_collect[j] = kappa;
    logpost = (p/2+alpha)*log(kappa)-kappa*(beta+dist);
    log_collect[j] = logpost;
  }

  Rcpp::NumericVector pctaccept = Rcpp::as<Rcpp::NumericVector>(numaccept) / iter;
  arma::vec Rr_best_match = best_match+1;
  arma::mat Rr_match_collect = match_collect+1;

  return (Rcpp::List::create(Rcpp::Named("best_match")=Rr_best_match,
                             Rcpp::Named("match_collect")=Rr_match_collect,
                             Rcpp::Named("dist_collect")=dist_collect,
                             Rcpp::Named("kappa_collect")=kappa_collect,
                             Rcpp::Named("log_collect")=log_collect,
                             Rcpp::Named("dist_min")=dist_min,
                             Rcpp::Named("pct_accept")=pctaccept));
}

Rcpp::List itercode(int iter, int n, int m, arma::vec mu_5,
                    arma::mat match_matrix, arma::mat qt_matrix,
                    arma::mat qt_fitted_matrix, int L, float tau, int times,
                    float kappa, float alpha, float beta, float powera,
                    arma::vec best_vec, arma::vec dist_vec,
                    arma::mat best_match_matrix, arma::vec mu_prior,
                    float var_const, arma::vec sumdist, int thin,
                    arma::mat mu_q, arma::mat mu_q_standard, float logmax,
                    int burnin, float AVG) {
  match_matrix = match_matrix - 1;
  best_match_matrix = best_match_matrix - 1;
  arma::rowvec scale(L);

  float increment,adjustcon,ratio,prob,u,n_dist,o_dist,dist,logpost,SigmaVar,rescale;
  int increment_int,newj,tempnew,tempold,tempx;
  arma::vec newmatchvec(3),oldmatchvec(3),idenmatchvec(3),qt2_5(m),match(L+1);
  arma::vec interynew(2*times),interyold(2*times),interx(2*times),scalevec(m),qt_5_fitted(m);
  Rcpp::IntegerVector Ixout(2*times),Ioriginal(L+1);
  arma::vec xout(2*times),internew(2*times),interold(2*times),mu_5_interx(2*times);
  arma::vec qt2_5_internew(2*times),qt2_5_interold(2*times),diff_ynew(2*times),diff_yold(2*times);
  arma::vec pnewvec(2),poldvec(2),original(L+1),idy(m);
  Rcpp::IntegerVector Irow = Rcpp::seq_len(m)-1;
  arma::vec row = Rcpp::as<arma::vec>(Irow);
  arma::vec karcher_res(m+1),tempkarcher_res(m+1),revscalevec(L*times),kappa_collect(iter);
  arma::vec log_collect(iter),MAP(m),Mean(m),qt2_5_warp(m),mu_5_warp(m);
  arma::mat tempmatch_matrix(L+1,n),bayes_warps(L+1,n);
  bayes_warps.fill(0);

  int q2LLlen = (m-1)*times+1;
  arma::mat q2LL_mat(q2LLlen,n);
  Rcpp::IntegerVector tempspan1(q2LLlen);
  tempspan1.fill(0);

  tempspan1 = Rcpp::seq_len(q2LLlen)-1;
  Rcpp::NumericVector temp_span2 = Rcpp::as<Rcpp::NumericVector>(tempspan1);
  float timesf = (float)(times);
  arma::vec q2LL_time = temp_span2*(1/timesf);
  Rcpp::IntegerVector q2L_time1 = Rcpp::seq_len(m)-1;
  Rcpp::NumericVector q2L_time2 = Rcpp::as<Rcpp::NumericVector>(q2L_time1);

  for (int LLL=0;LLL<n;LLL++ )
  {
    q2LL_mat.col(LLL) = approx(m,q2L_time2,qt_matrix.col(LLL),q2LLlen,q2LL_time);
  }


  Rcpp::RNGScope scope;

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
          increment =  Rcpp::as<float>(Rcpp::rnorm(1,0,tau));
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
            Ixout = Rcpp::seq_len(2*times)+times*(i-1)-1;
            xout = Rcpp::as<arma::vec>(Ixout);
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
            u =  Rcpp::as<float>(Rcpp::runif(1));
            match[i] = (u < prob)?(newj):(match[i]);
          }
        }
      }
      Ioriginal = (Rcpp::seq_len(L+1)-1)*times;
      original = Rcpp::as<arma::vec>(Ioriginal);
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
      mu_q.row(j/thin) = arma::conv_to<arma::rowvec>::from(mu_5);
      karcher_res = findinv(match_matrix,times);
      for (int jj=0;jj<m;jj++){revscalevec[jj] = sqrt(karcher_res(jj+1)-karcher_res(jj));}
      mu_5_warp =  approx(m,q2L_time2,mu_5,m,karcher_res);
      for (int qq=0;qq<m;qq++){mu_q_standard(j/thin,qq) = revscalevec[qq]*mu_5_warp[qq];}
    }
    kappa = Rcpp::as<float>(Rcpp::rgamma(1,n*m/2+alpha, 1/(sum(dist_vec)+beta)));
    kappa_collect[j] = kappa;
    logpost = (n*m/2+alpha)*log(kappa)-kappa*(beta+sum(dist_vec));
    log_collect[j] = logpost;
    if(logpost > logmax){logmax = logpost;MAP = mu_5;}
    Mean = (var_const*m)/(m+2*kappa*n*var_const)*(2*kappa*n*mean(qt_fitted_matrix,1)/m +mu_prior/var_const);
    SigmaVar = (var_const*m)/(m+2*kappa*n*var_const);
    mu_5 = Mean + Rcpp::as<arma::vec>(Rcpp::rnorm(m,0,sqrt(SigmaVar)));
    rescale = sqrt(m/sum(square(mu_5)));
    mu_5 = rescale*mu_5;
  }

  arma::mat Rr_best_match_matrix = best_match_matrix+1;
  arma::mat R_bayes_warps = bayes_warps+1;

  return Rcpp::List::create(Rcpp::Named("mu.q.standard") = mu_q_standard,
                            Rcpp::Named("mu.q") = mu_q,
                            Rcpp::Named("best_match.matrix") = Rr_best_match_matrix,
                            Rcpp::Named("kappafamily") = kappa_collect,
                            Rcpp::Named("sumdist") = sumdist,
                            Rcpp::Named("log.posterior") = log_collect,
                            Rcpp::Named("dist") = sum(best_vec),
                            Rcpp::Named("MAP") = MAP,
                            Rcpp::Named("bayes_warps") = R_bayes_warps);
}

arma::vec approx(int nd, arma::vec xd, arma::vec yd, int ni, arma::vec xi)
{
  int i,k;
  double t;
  arma::vec yi(ni);
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

arma::vec findinv(arma::mat warps, int times)
{

  int m = warps.n_rows;
  int n = warps.n_cols;
  arma::mat psi_m(m-1,n);
  psi_m.zeros();

  for(int j = 0; j<n; j++)
  {
    for (int hh = 0; hh< m-1 ; hh++)
    {
      psi_m(hh,j) = sqrt( (warps(hh+1,j)-warps(hh,j))/times);
    }
  }

  arma::vec w = mean(psi_m,1);
  arma::vec mupsi = w/arma::as_scalar(sqrt(sum(pow(w,2)/(m-1))));
  arma::mat v_m(m-1,n);
  v_m.zeros();
  arma::vec mupsi_update(m-1);
  float check = 1;

  while (check > 0.01)
  {
    for (int i = 0;i<n;i++)
    {
      double theta = acos(sum(mupsi.t()*psi_m.col(i)/(m-1)));
      v_m.col(i) = theta/sin(theta)*(psi_m.col(i)-cos(theta)*mupsi);
    }
    arma::vec vbar = mean(v_m,1);
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
  arma::vec karcher_s = cumsum(pow(mupsi_update,2))*times;
  karcher_s.insert_rows(0,1);
  Rcpp::IntegerVector Irow = (Rcpp::seq_len(m)-1)*times;
  arma::vec row = Rcpp::as<arma::vec>(Irow);
  Rcpp::IntegerVector Iout = Rcpp::seq_len((m-1)*times+1)-1;
  arma::vec out = Rcpp::as<arma::vec>(Iout);
  arma::vec invidy = approx(m,karcher_s,row,(m-1)*times+1,out);

  return invidy;

}

arma::vec R_diff(arma::vec x)
{
  int n = x.n_elem;
  arma::vec x_diff(n-1);

  for (int i=0;i<(n-1);i++)
  {
    x_diff[i] = x[i+1]-x[i];
  }

  return x_diff;

}
