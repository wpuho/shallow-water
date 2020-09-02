function phiuv = timeintegral_h(h_H,hM,u,v,nx,my,dx,dy,dt)

%%
phiu = zeros(my,nx);
h_hM = h_H - hM;

for j = 1:my
    for i = 2:nx-1

        if (u(j,i) >= 0)      
             if (i==2 || i==nx-1)
                 phiu(j,i) = h_H(j,i) - dt*u(j,i)*((h_hM(j,i)-h_hM(j,i-1))/(dx));
             else
                 phiu(j,i) = h_H(j,i) - dt*u(j,i)*((h_hM(j,i+1)-h_hM(j,i-1))/(2*dx) - (h_hM(j,i+1)-3*h_hM(j,i)+3*h_hM(j,i-1)-h_hM(j,i-2))/(6*dx));
             end
             
        elseif (u(j,i) < 0)
             if (i==nx-1 || i==2)
                 phiu(j,i) = h_H(j,i) - dt*u(j,i)*((h_hM(j,i+1)-h_hM(j,i))/(dx));
             else
                 phiu(j,i) = h_H(j,i) - dt*u(j,i)*((h_hM(j,i+1)-h_hM(j,i-1))/(2*dx) - (h_hM(j,i+2)-3*h_hM(j,i+1)+3*h_hM(j,i)-h_hM(j,i-1))/(6*dx));
             end
        end
        
    end
end

%
phiu(:,1) = phiu(:,2);
phiu(:,nx) = phiu(:,nx-1);
%%
phiuv = zeros(my,nx);
h_hM = phiu - hM;

for i = 1:nx
    for j = 2:my-1

            if (v(j,i) >= 0)
                if (j==2 || j==my-1)
                    phiuv(j,i) = phiu(j,i) - dt*v(j,i)*((h_hM(j,i)-h_hM(j-1,i))/(dy));
                else
                    phiuv(j,i) = phiu(j,i) - dt*v(j,i)*((h_hM(j+1,i)-h_hM(j-1,i))/(2*dy) - (h_hM(j+1,i)-3*h_hM(j,i)+3*h_hM(j-1,i)-h_hM(j-2,i))/(6*dy));
                end
            elseif (v(j,i) < 0)
                if (j==my-1 || j==2)
                    phiuv(j,i) = phiu(j,i) - dt*v(j,i)*((h_hM(j+1,i)-h_hM(j,i))/(dy));
                else
                    phiuv(j,i) = phiu(j,i) - dt*v(j,i)*((h_hM(j+1,i)-h_hM(j-1,i))/(2*dy) - (h_hM(j+2,i)-3*h_hM(j+1,i)+3*h_hM(j,i)-h_hM(j-1,i))/(6*dy));
                end
            end        
        
    end
end
%
phiuv(1,:) = phiuv(2,:);
phiuv(my,:) = phiuv(my-1,:);

phiuv(1,1) = (phiuv(2,1)+phiuv(1,2))/2; %phi(2,2);
phiuv(1,nx) = (phiuv(1,nx-1)+phiuv(2,nx))/2;%phi(2,nx-1);
phiuv(my,1) = (phiuv(my-1,1)+phiuv(my,2))/2;%phi(my-1,2);
phiuv(my,nx) = (phiuv(my,nx-1)+phiuv(my-1,nx))/2;%phi(my-1,nx-1);
%phiuv = zerogradientBC2(phiuv,nx,my);
