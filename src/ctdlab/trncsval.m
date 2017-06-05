function A = trncsval(A,tol)

    ikeep = find(abs(A.lambda)>tol);
    
    if isempty(ikeep)
      disp('trncsval.m:  all svals are below tolerance, not truncating\n')
      return
    end
    
    A = getterms(A,ikeep);

end