function obj =  val_test (tau) 
      tic;
      dlmwrite('tau.txt', tau);
      system('/usr/local/bin/R CMD BATCH test_2.R');
      load aveY.txt
      obj = aveY;
      disp(obj)
      toc;
end
    
 