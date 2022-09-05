% semplifica i valori di C(i,i)
%del tipo res=simplify_this_C(temp,dq,N)
function res=simplify_this_C(temp,dq,N)
            tempsimp=0;
            for k=1:N
                for h=1:N
                    debug=[k,h];
                    if (k<h)
                        tempsimp=tempsimp+simplify(expand(diff(diff(temp,dq(k)),dq(h))))*dq(k)*dq(h);
                    end
                    if (k==h)
                        tempsimp=tempsimp+simplify(expand(diff(diff(temp,dq(k)),dq(h))))/2*dq(k)*dq(h);
                    end
                end
            end
            res=tempsimp;
end
