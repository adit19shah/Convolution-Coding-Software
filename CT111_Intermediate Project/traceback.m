% Traceback starting from s_1 i.e. [0 0] present in the last column
function received=traceback(s,K,k,path_m,branch_m)
state=s(1);
% Initializing the vector for final decoded sequence with '-1's
received=-ones(1,k+K-1);    

if(K==3)
    for b=1:k+2
      switch state
        case s(1)
          received(1,k+3-b)=0; 
          % tracing back towards the branch which offers minimum sum of
          % path metric and branch metric
          if(min(path_m(3,k+3-b)+branch_m(6,k+3-b),path_m(4,k+3-b)+branch_m(8,k+3-b))==path_m(3,k+3-b)+branch_m(6,k+3-b))
            state=s(2);
          else
            state=s(1);
          end

        case s(2)
          received(1,k+3-b)=0;
          if(min(path_m(1,k+3-b)+branch_m(2,k+3-b),path_m(2,k+3-b)+branch_m(4,k+3-b))==path_m(1,k+3-b)+branch_m(2,k+3-b))
            state=s(4);
          else
            state=s(3);
          end

         case s(3)
          received(1,k+3-b)=1;
          if(min(path_m(3,k+3-b)+branch_m(5,k+3-b),path_m(4,k+3-b)+branch_m(7,k+3-b))==path_m(3,k+3-b)+branch_m(5,k+3-b))
            state=s(2);
          else
            state=s(1);
          end

         case s(4)
          received(1,k+3-b)=1;
          if(min(path_m(1,k+3-b)+branch_m(1,k+3-b),path_m(2,k+3-b)+branch_m(3,k+3-b))==path_m(1,k+3-b)+branch_m(1,k+3-b))
            state=s(4);
          else
            state=s(3);
          end   
      end
    end

elseif(K==4)
    for b=1:k+3
      switch state
        case s(1)
          received(1,k+4-b)=0;
          % tracing the path on the basis of minimum value of path metric connected
          % to this state
          if(min(path_m(7,k+4-b)+branch_m(14,k+4-b),path_m(8,k+4-b)+branch_m(16,k+4-b))==path_m(7,k+4-b)+branch_m(15,k+4-b))
            state=s(2);
          else
            state=s(1);
          end

        case s(2)
          received(1,k+4-b)=0;
          if(min(path_m(5,k+4-b)+branch_m(10,k+4-b),path_m(6,k+4-b)+branch_m(12,k+4-b))==path_m(5,k+4-b)+branch_m(10,k+4-b))
            state=s(4);
          else
            state=s(3);
          end

         case s(3)
          received(1,k+4-b)=0;
          if(min(path_m(3,k+4-b)+branch_m(6,k+4-b),path_m(4,k+4-b)+branch_m(8,k+4-b))==path_m(3,k+4-b)+branch_m(6,k+4-b))
            state=s(6);
          else
            state=s(5);
          end

         case s(4)
          received(1,k+4-b)=0;
          if(min(path_m(1,k+4-b)+branch_m(2,k+4-b),path_m(2,k+4-b)+branch_m(4,k+4-b))==path_m(1,k+4-b)+branch_m(2,k+4-b))
            state=s(8);
          else
            state=s(7);
          end 

         case s(5)
          received(1,k+4-b)=1;
          if(min(path_m(7,k+4-b)+branch_m(13,k+4-b),path_m(8,k+4-b)+branch_m(15,k+4-b))==path_m(7,k+4-b)+branch_m(13,k+4-b))
            state=s(2);
          else
            state=s(1);
          end

         case s(6)
          received(1,k+4-b)=1;
          if(min(path_m(5,k+4-b)+branch_m(9,k+4-b),path_m(6,k+4-b)+branch_m(11,k+4-b))==path_m(5,k+4-b)+branch_m(9,k+4-b))
            state=s(4);
          else
            state=s(3);
          end

         case s(7)
          received(1,k+4-b)=1;
          if(min(path_m(3,k+4-b)+branch_m(5,k+4-b),path_m(4,k+4-b)+branch_m(7,k+4-b))==path_m(3,k+4-b)+branch_m(5,k+4-b))
            state=s(6);
          else
            state=s(5);
          end

         case s(8)
          received(1,k+4-b)=1;
          if(min(path_m(1,k+4-b)+branch_m(1,k+4-b),path_m(2,k+4-b)+branch_m(3,k+4-b))==path_m(1,k+4-b)+branch_m(1,k+4-b))
            state=s(8);
          else
            state=s(7);
          end
      end
    end
end


