% ------------------------------------------------------------------------%

% filename: Ternary_1DFD_gen.m

% This function generates Ternary piecewise-constant signal using the Gibb's

% sampling method. 

% JWKANG 2013, SEP. (jwkkang@gist.ac.kr), updated at 2014 Nov

%-------------------------------------------------------------------------%

function X=Ternary_1DFD_gen(N,q,sigma0)


    while(1)
        Numofsupp=length(find(rand(1,N)<q));
        Xstate=randerr(1,N,Numofsupp);
        supp=find(Xstate==1);
        S=zeros(N,1);
        S(supp)=1;
        Xdiff = S.*(floor(3*rand(N,1))-1)*sigma0;
        X = cumsum(Xdiff);
        
        if sum(diff(X))~=0
            break;
        end
    end

    




end

