function [ pry ] = estPrY(y, ml, cl, rl,  sd, n) 
   pry = 0.00001;
   for  i = 1:n
       m = ml( i );
       c = cl( i );
       r = rl( i );
       t = ( r > 0 ) - ( r <= 0 );
      % a = al( i );
      % calculate   prY = cdfY( y, m, c, r)
      % create a list of each y  with the same value and same dimension of
      % m c r 
      % I did not include variance here
      % pry = pry +  normcdf( (y - m - t * c)/sd );
%      pry = pry +  normcdf( (y - m - t * c)/sd );
       pry = pry +  1/2 * (1 + erf ( (y - m - t * c)/ (sqrt(2)*sd) ) ) ;
   end
   pry = pry / n ;
end