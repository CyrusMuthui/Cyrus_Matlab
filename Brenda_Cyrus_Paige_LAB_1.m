%{
Brenda Cheptoo - ENE211-0004/2018
Cyrus Muthui   - ENE211-0010/2018
Paige Muva     - ENE211-0017/2018
%}
i=int64(0);
V1= 1;
delta_1=0;
V2=1;
delta_2= 0;
V3=1;
delta_3= 0;
V4=1;
delta_4= 0;
V5=1.02;
delta_5= 0;

p2=-0.96;
q2=-0.62;
p3=-0.35;
q3=-0.14;
p4=-0.16;
q4=-0.08;
p5=0.24;
q5=-0.11;


Y=[2.6923-13.4096i  -1.9231+9.6154i   0    0   -0.7692+3.8462i ;-1.9231+9.6154i   3.6541-18.1924i  -0.962+4.808i   0   -0.769+3.846i; 0   -0.962+4.808i   2.2115-11.0023i   -0.7692+3.846i   -0.4808+2.4038i; 0   0  -0.7692+3.846i  1.1538-5.6741i   -0.3846+1.9231i;-0.7692+3.8462i  -0.769+3.846i  -0.4808+2.438i  -0.3846+1.9231i  2.4038-11.8942i];
jac=[];
del_f=[];
del_x=[];

fprintf('      Iter        V2     V3       V4      Delta_2     Delta_3      Delta_4     Delta_5     Q5    \n\n' );

for i= 0:1:3
    
p2_calc= V1*V2*abs(Y(1,2))*cos(angle(Y(1,2))-delta_2 +delta_1)+ power(V2,2)*abs(Y(2,2))*cos(angle(Y(2,2)))+ V3*V2*abs(Y(2,3))*cos(angle(Y(2,3))-delta_2 +delta_3)+ V4*V2*abs(Y(2,4))*cos(angle(Y(2,4))-delta_2 +delta_4)+V5*V2*abs(Y(2,5))*cos(angle(Y(2,5))-delta_2 +delta_5);
q2_calc= -V1*V2*abs(Y(1,2))*sin(angle(Y(1,2))-delta_2 +delta_1)- power(V2,2)*abs(Y(2,2))*sin(angle(Y(2,2)))-V3*V2*abs(Y(2,3))*sin(angle(Y(2,3))-delta_2 +delta_3)- V4*V2*abs(Y(2,4))*sin(angle(Y(2,4))-delta_2 +delta_4)-V5*V2*abs(Y(2,5))*sin(angle(Y(2,5))-delta_2 +delta_5);
p3_calc= V3*V1*abs(Y(3,1))*cos(angle(Y(3,1))-delta_3 +delta_1)+ power(V3,2)*abs(Y(3,3))*cos(angle(Y(3,3)))+ V3*V2*abs(Y(3,2))*cos(angle(Y(3,2))-delta_3 +delta_2)+ V3*V4*abs(Y(3,4))*cos(angle(Y(3,4))-delta_3 +delta_4)+V3*V5*abs(Y(3,5))*cos(angle(Y(3,5))-delta_3 +delta_5);
q3_calc= -V3*V1*abs(Y(3,1))*sin(angle(Y(3,1))-delta_3 +delta_1)- power(V3,2)*abs(Y(3,3))*sin(angle(Y(3,3)))-V3*V2*abs(Y(3,2))*sin(angle(Y(3,2))-delta_3 +delta_2)- V3*V4*abs(Y(3,4))*sin(angle(Y(3,4))-delta_3 +delta_4)-V3*V5*abs(Y(3,5))*sin(angle(Y(3,5))-delta_3 +delta_5);
p4_calc= V4*V1*abs(Y(4,1))*cos(angle(Y(4,1))-delta_4 +delta_1)+ power(V4,2)*abs(Y(4,4))*cos(angle(Y(4,4)))+ V4*V2*abs(Y(4,2))*cos(angle(Y(4,2))-delta_4 +delta_2)+ V4*V3*abs(Y(4,3))*cos(angle(Y(4,3))-delta_4 +delta_3)+V4*V5*abs(Y(4,5))*cos(angle(Y(4,5))-delta_4 +delta_5);
q4_calc= -V4*V1*abs(Y(4,1))*sin(angle(Y(4,1))-delta_4 +delta_1)- power(V4,2)*abs(Y(4,4))*sin(angle(Y(4,4)))- V4*V2*abs(Y(4,2))*sin(angle(Y(4,2))-delta_4 +delta_2)- V4*V3*abs(Y(4,3))*sin(angle(Y(4,3))-delta_4 +delta_3)-V4*V5*abs(Y(4,5))*sin(angle(Y(4,5))-delta_4 +delta_5);
p5_calc=  V5*V1*abs(Y(5,1))*cos(angle(Y(5,1))-delta_5 +delta_1)+ power(V5,2)*abs(Y(5,5))*cos(angle(Y(5,5)))+ V5*V2*abs(Y(5,2))*cos(angle(Y(5,2))-delta_5 +delta_2)+ V5*V3*abs(Y(5,3))*cos(angle(Y(5,3))-delta_5 +delta_3)+V5*V4*abs(Y(5,4))*cos(angle(Y(5,4))-delta_5 +delta_4);


                 
    jac(1, 1)= V1*V2*abs(Y(2,1))*sin(angle(Y(1,2))-delta_2 +delta_1)+V3*V2*abs(Y(2,3))*sin(angle(Y(2,3))-delta_2 +delta_3)+ V4*V2*abs(Y(2,4))*sin(angle(Y(2,4))-delta_2 +delta_4)+V5*V2*abs(Y(2,5))*sin(angle(Y(2,5))-delta_2 +delta_5);
    jac(1, 2)= - V3*V2*abs(Y(2,3))*sin(angle(Y(2,3))-delta_2 +delta_3);
    jac(1, 3)=  -V4*V2*abs(Y(2,4))*sin(angle(Y(2,4))-delta_2 +delta_4); 
    jac(1, 4)= - V5*V2*abs(Y(2,5))*sin(angle(Y(2,5))-delta_2 +delta_5);
    jac(1, 5)= V1*abs(Y(1,2))*cos(angle(Y(1,2))-delta_2 +delta_1)+ 2*V2*abs(Y(2,2))*cos(angle(Y(2,2)))+ V3*abs(Y(2,3))*cos(angle(Y(2,3))-delta_2 +delta_3)+ V4*abs(Y(2,4))*cos(angle(Y(2,4))-delta_2 +delta_4)+V5*abs(Y(2,5))*cos(angle(Y(2,5))-delta_2 +delta_5);
    jac(1, 6)=  V2*abs(Y(2,3))*cos(angle(Y(2,3))-delta_2 +delta_3);
    jac(1, 7)=  V2*abs(Y(2,4))*cos(angle(Y(2,4))-delta_2 +delta_4);
    jac(2, 1)= - V3*V2*abs(Y(3,2))*sin(angle(Y(3,2))-delta_3 +delta_2);
    jac(2, 2)=  V3*V1*abs(Y(3,1))*sin(angle(Y(3,1))-delta_3 +delta_1)+V3*V2*abs(Y(3,2))*sin(angle(Y(3,2))-delta_3 +delta_2)+ V3*V4*abs(Y(3,4))*sin(angle(Y(3,4))-delta_3 +delta_4)+V3*V5*abs(Y(3,5))*sin(angle(Y(3,5))-delta_3 +delta_5);
    jac(2, 3)=  -V3*V4*abs(Y(3,4))*sin(angle(Y(3,4))-delta_3 +delta_4);
    jac(2, 4)=  -V3*V5*abs(Y(3,5))*sin(angle(Y(3,5))-delta_3 +delta_5);
    jac(2, 5)=   V3*abs(Y(3,2))*cos(angle(Y(3,2))-delta_3 +delta_2);
    jac(2, 6)=  V3*abs(Y(3,1))*cos(angle(Y(3,1))-delta_3 +delta_1)+ 2*V3*abs(Y(3,3))*cos(angle(Y(3,3)))+ V3*abs(Y(3,2))*cos(angle(Y(3,2))-delta_3 +delta_2)+V4*abs(Y(3,4))*cos(angle(Y(3,4))-delta_3 +delta_4)+V5*abs(Y(3,5))*cos(angle(Y(3,5))-delta_3 +delta_5);
    jac(2, 7)=   V3*abs(Y(3,4))*cos(angle(Y(3,4))-delta_3 +delta_4);
    jac(3, 1)=  -V4*V2*abs(Y(4,2))*sin(angle(Y(4,2))-delta_4 +delta_2);
    jac(3, 2)=  - V4*V3*abs(Y(4,3))*sin(angle(Y(4,3))-delta_4 +delta_3);
    jac(3, 3)=  V4*V1*abs(Y(4,1))*sin(angle(Y(4,1))-delta_4 +delta_1)+ V4*V2*abs(Y(4,2))*sin(angle(Y(4,2))-delta_4 +delta_2)+ V4*V3*abs(Y(4,3))*sin(angle(Y(4,3))-delta_4 +delta_3)+V4*V5*abs(Y(4,5))*sin(angle(Y(4,5))-delta_4 +delta_5);
    jac(3, 4)=  -V4*V5*abs(Y(4,5))*sin(angle(Y(4,5))-delta_4 +delta_5);
    jac(3, 5)=  V4*abs(Y(4,2))*cos(angle(Y(4,2))-delta_4 +delta_2);
    jac(3, 6)= V4*abs(Y(4,3))*cos(angle(Y(4,3))-delta_4 +delta_3);
    jac(3, 7)=  V1*abs(Y(4,1))*cos(angle(Y(4,1))-delta_4 +delta_1)+ 2*V4*abs(Y(4,4))*cos(angle(Y(4,4)))+ V2*abs(Y(4,2))*cos(angle(Y(4,2))-delta_4 +delta_2)+V3*abs(Y(4,3))*cos(angle(Y(4,3))-delta_4 +delta_3)+V5*abs(Y(4,5))*cos(angle(Y(4,5))-delta_4 +delta_5);
    jac(4, 1)= -V5*V2*abs(Y(5,2))*sin(angle(Y(5,2))-delta_5 +delta_2);
    jac(4, 2)=  -V5*V3*abs(Y(5,3))*sin(angle(Y(5,3))-delta_5 +delta_3);
    jac(4, 3)=  -V5*V4*abs(Y(5,4))*sin(angle(Y(5,4))-delta_5 +delta_4);
    jac(4, 4)=  V5*V1*abs(Y(5,1))*sin(angle(Y(5,1))-delta_5 +delta_1)+ V5*V2*abs(Y(5,2))*sin(angle(Y(5,2))-delta_5 +delta_2)+ V5*V3*abs(Y(5,3))*sin(angle(Y(5,3))-delta_5 +delta_3)+V5*V4*abs(Y(5,4))*sin(angle(Y(5,4))-delta_5 +delta_4);
    jac(4, 5)=  V5*abs(Y(5,2))*cos(angle(Y(5,2))-delta_5 +delta_2);
    jac(4, 6)= V5*abs(Y(5,3))*cos(angle(Y(5,3))-delta_5 +delta_3);
    jac(4, 7)= V5*abs(Y(5,4))*cos(angle(Y(5,4))-delta_5 +delta_4);
    jac(5, 1)=  V1*V2*abs(Y(1,2))*cos(angle(Y(1,2))-delta_2 +delta_1)+V3*V2*abs(Y(2,3))*cos(angle(Y(2,3))-delta_2 +delta_3)+ V4*V2*abs(Y(2,4))*cos(angle(Y(2,4))-delta_2 +delta_4)+V5*V2*abs(Y(2,5))*cos(angle(Y(2,5))-delta_2 +delta_5);
    jac(5, 2)= -V3*V2*abs(Y(2,3))*cos(angle(Y(2,3))-delta_2 +delta_3);
    jac(5, 3)=  -V4*V2*abs(Y(2,4))*cos(angle(Y(2,4))-delta_2 +delta_4);
    jac(5, 4)= -V5*V2*abs(Y(2,5))*cos(angle(Y(2,5))-delta_2 +delta_5);
    jac(5, 5)=  -V1*abs(Y(1,2))*sin(angle(Y(1,2))-delta_2 +delta_1)- 2*V2*abs(Y(2,2))*sin(angle(Y(2,2)))-V3*abs(Y(2,3))*sin(angle(Y(2,3))-delta_2 +delta_3)- V4*abs(Y(2,4))*sin(angle(Y(2,4))-delta_2 +delta_4)-V5*abs(Y(2,5))*sin(angle(Y(2,5))-delta_2 +delta_5);
    jac(5, 6)= -V2*abs(Y(2,3))*sin(angle(Y(2,3))-delta_2 +delta_3);
    jac(5, 7)=- V2*abs(Y(2,4))*sin(angle(Y(2,4))-delta_2 +delta_4);
    jac(6, 1)= - V3*V2*abs(Y(3,2))*cos(angle(Y(3,2))-delta_3 +delta_2);
    jac(6, 2)= V3*V1*abs(Y(3,1))*cos(angle(Y(3,1))-delta_3 +delta_1)+V3*V2*abs(Y(3,2))*cos(angle(Y(3,2))-delta_3 +delta_2)+ V3*V4*abs(Y(3,4))*cos(angle(Y(3,4))-delta_3 +delta_4)+V3*V5*abs(Y(3,5))*cos(angle(Y(3,5))-delta_3 +delta_5);
    jac(6, 3)= -V3*V4*abs(Y(3,4))*cos(angle(Y(3,4))-delta_3 +delta_4);
    jac(6, 4)= -V3*V5*abs(Y(3,5))*cos(angle(Y(3,5))-delta_3 +delta_5);
    jac(6, 5)= -V3*abs(Y(3,2))*sin(angle(Y(3,2))-delta_3 +delta_2);
    jac(6, 6)= -V1*abs(Y(3,1))*sin(angle(Y(3,1))-delta_3 +delta_1)- 2*V3*abs(Y(3,3))*sin(angle(Y(3,3)))-V2*abs(Y(3,2))*sin(angle(Y(3,2))-delta_3 +delta_2)- V3*abs(Y(3,4))*sin(angle(Y(3,4))-delta_3 +delta_4)-V5*abs(Y(3,5))*sin(angle(Y(3,5))-delta_3 +delta_5);
    jac(6, 7)= - V3*abs(Y(3,4))*sin(angle(Y(3,4))-delta_3 +delta_4);
    jac(7, 1)= - V4*V2*abs(Y(4,2))*cos(angle(Y(4,2))-delta_4 +delta_2);
    jac(7, 2)=  - V4*V3*abs(Y(4,3))*cos(angle(Y(4,3))-delta_4 +delta_3);
    jac(7, 3)=  V4*V1*abs(Y(4,1))*cos(angle(Y(4,1))-delta_4 +delta_1)+ V4*V2*abs(Y(4,2))*cos(angle(Y(4,2))-delta_4 +delta_2)+ V4*V3*abs(Y(4,3))*cos(angle(Y(4,3))-delta_4 +delta_3)+V4*V5*abs(Y(4,5))*cos(angle(Y(4,5))-delta_4 +delta_5);
    jac(7, 4)=  -V4*V5*abs(Y(4,5))*cos(angle(Y(4,5))-delta_4 +delta_5);
    jac(7, 5)= - V4*abs(Y(4,2))*sin(angle(Y(4,2))-delta_4 +delta_2);
    jac(7, 6)= -V4*abs(Y(4,3))*sin(angle(Y(4,3))-delta_4 +delta_3);
    jac(7, 7)= -V1*abs(Y(4,1))*sin(angle(Y(4,1))-delta_4 +delta_1)- 2*V4*abs(Y(4,4))*sin(angle(Y(4,4)))- V2*abs(Y(4,2))*sin(angle(Y(4,2))-delta_4 +delta_2)- V3*abs(Y(4,3))*sin(angle(Y(4,3))-delta_4 +delta_3)-V5*abs(Y(4,5))*sin(angle(Y(4,5))-delta_4 +delta_5);
    
              
    del_f(1,1)= p2-p2_calc;
     del_f(2,1)= p3-p3_calc;
     del_f(3,1)= p4-p4_calc;
      del_f(4,1)= p5-p5_calc;
      del_f(5,1)= q2-q2_calc;
     del_f(6,1)= q3-q3_calc;
     del_f(7,1)= q4-q4_calc;
    
    
    del_x= inv(jac)*del_f;
    
    delta_2= delta_2 +del_x(1,1);
    delta_3= delta_3 +del_x(2,1);
    delta_4= delta_4 +del_x(3,1);
    delta_5= delta_5 +del_x(4,1);
      V2= V2+ del_x(5,1);
      V3= V3+ del_x(6,1);
      V4= V4+ del_x(7,1);
    
    disp([i,  V2 ,   V3,  V4  ,  rad2deg(delta_2), rad2deg(delta_3), rad2deg(delta_4), rad2deg(delta_5)  q5 ]);
end