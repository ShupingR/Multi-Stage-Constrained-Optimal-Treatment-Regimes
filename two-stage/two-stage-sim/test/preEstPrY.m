% ? how to deal with d(x1).
% usued variable dont define.
% H2 depend on d1 as well.
% normalize d1 and d2
function [ meanYEst ] = preEstPrY (tau, Y, H2, A2, H1, A1, n, sign)
  %================================================
  % @Fun: to estimation marginal cumulative distribution of Y under any  
  % treatment regime d = (d_1, d_2) in the class of linear decision rule    
  % @Input:                                                                                               
  %     tau: linear decision rule indices to be max over                                     
  %     Y H2 A2 H1 A1: data for dtr                                                          
  %     n : resapme size               
  %     sign: the sign of the value to be returned, choose from -1 or 1
  %================================================
 
    tau1 = tau( 1 : 2 );
    %/norm(tau( 1 : 2 ));
    tau2 = tau( 3 : 4 );
    %/norm(tau( 3 : 4 ));
    %%
    %-------------------------------
    %  Estimate all the parameters  
    %-------------------------------
    
    %% Estimate mY cY
    xdata2 = [ H2 , A2 .* H2( :, 1), A2 .* H2( :, 2) ];
    hatYbeta= xdata2 \ Y;
    hatYbeta0 = hatYbeta( 1 : 2 ); %*
    hatYbeta1 = hatYbeta( 3 : 4 ); %*
    hatMY = H2 * hatYbeta0;
    hatCY = H2 * hatYbeta1;
    
    %% Estimated std of residuals of eY
    hatErrY = Y - hatMY - A2 .* hatCY;
    hatsigmaEY = std( hatErrY ); %* 
        % homogeneous assumed 
        % ErrY/stdErrY ~ Normal (0, 1) assumed

    % first stage dataset to be regressed on for estimating parameters
    xdata1 = [ H1, A1 .* H1(:, 1), A1 .* H1(:, 2) ]; 
 
    %% residuals of mY (mean and variance modelling)
    hatmuMYgamma = xdata1 \ hatMY ; %*
    hatmuMY = xdata1 * hatmuMYgamma;
    hatsigmaMYgamma = xdata1 \ ( 2 * log ( abs( hatMY - hatmuMY ) ) ) ; %*
    hatsigmaMY =  exp( xdata1 * hatsigmaMYgamma / 2 ) .^ 0.5 ; 
    hatErrStdMY = ( hatMY - hatmuMY ) ./ hatsigmaMY ;%*
    
    %% residuals of cY
    hatmuCYgamma = xdata1 \ hatCY ; %*
    hatmuCY = xdata1 * hatmuCYgamma ;
    hatsigmaCYgamma = xdata1 \ ( 2 * log ( abs ( hatCY - hatmuCY ) ) ) ;  %*
    hatsigmaCY =  exp( xdata1 * hatsigmaCYgamma / 2 ) .^ 0.5 ; 
    hatErrStdCY = ( hatCY - hatmuCY ) ./ hatsigmaCY ; %*
    
    %% residuals of r2
    r2 = H2 * tau2' ;
    hatmuR2gamma = xdata1 \ r2 ; %*
    hatmuR2 = xdata1 * hatmuR2gamma ;
    hatsigmaR2gamma = xdata1 \ ( 2 * log ( abs ( r2 - hatmuR2 ) ) ) ; %*
    hatsigmaR2 =  exp( xdata1 * hatsigmaR2gamma / 2 ) .^ 0.5 ; 
    hatErrStdR2 = ( r2 - hatmuR2 ) ./ hatsigmaR2 ;
    
    %%
    %----------------------------------------
    %  Sample eY, mY, cY under a regime d 
    %  using the estimated parameters         
    %----------------------------------------
    % xdatad2 = [ H2 , d2 .* H2( :, 1), d2 .* H2( :, 2) ];
    
    %% estimated first stage data and its transform that follows d1
    r1 = H1 * tau1';
    d1 = ( r1 > 0 ) - ( r1 <= 0 );
    xdatad1 = [ H1, d1 .* H1(:, 1), d1 .* H1(:, 2) ];
    
    hatmuMYd1 = xdatad1 * hatmuMYgamma;
    hatsigmaMYd1 =  exp( xdatad1 * hatsigmaMYgamma / 2 ) .^ 0.5 ; 
    
    hatmuCYd1 = xdatad1 * hatmuCYgamma;
    hatsigmaCYd1 =  exp( xdatad1 * hatsigmaCYgamma / 2 ) .^ 0.5 ; 
    
    % hatmuR2d1 = hatmuR2;
    % hatsigmaR2 = hatsigmaR2;
    
%   rng(1);
    %% get a sample of m c r under regime d
    rng(1,'twister');
    rErr = copula( hatErrStdMY, hatErrStdCY, hatErrStdR2, n) ;
        % use copula to resample from the standardized residuals of mY, cY,
        % R2, as all the dependence is catched in mean variance modeling of
        % muMY, sigmaMY, muCY and sigmaCY. The std residuals is distributed
        % the same under any frist stage rule d1 and with any first-stage covariates H1
    mRs = rErr( : , 1) .*  hatsigmaMYd1 + hatmuMYd1 ;
    cRs= rErr( : , 2) .*  hatsigmaCYd1  + hatmuCYd1 ;
    rRs = rErr( : , 3) .*  hatsigmaR2  + hatmuR2 ;
        % no need to estimated H2(H1, d1), 
        % as we estimated and are intergating wrt m(H2) | H1, d1
        
    %%
    prYlist  = zeros(50, 1);
    yleft = linspace( ceil(min(Y)) , floor(max(Y)) , 51 );
    yright = linspace( ceil(min(Y)) , floor(max(Y)) , 51 );
    ymid=  ( yleft( 1 : end-1 )+ yright( 2 : end ) ) / 2;
    %create a list of Y ranges
    for i  = 1: 50
        % prYlist(i)  = estPrY( yright( i + 1 ), mRs, cRs, rRs, hatsigmaEY, n ) - estPrY( yleft( i ), mRs, cRs, rRs, hatsigmaEY, n );
        prYlist(i)  = estPrY( yright( i + 1 ), mRs, cRs, rRs,  hatsigmaEY, n ) - estPrY( yleft( i ), mRs, cRs, rRs,  hatsigmaEY, n );
    end


    % generate from the estimated prYlist using multinomial sampling, sample
    % size M 
    m = 10000;
    p = prYlist / sum(prYlist) ;
    rng(2,'twister');
    ynum = mnrnd(m, p, 1);
    % project it back to Y scale
    % yRs = zeros ( length(p), 2 );
    %Est = [];
    %for i = 1: length(p)
    %    yEst = [ yEst, repmat( ymid( i ) , 1, ynum( i ) ) ];
    %end
    meanYEst =  sign *  sum( ymid .* ynum ) / sum (ynum);
end