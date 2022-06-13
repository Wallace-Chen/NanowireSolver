% piezoeletric polarization longitudinal in a wurtzite crystal
% theta is the angle between the natural c-axis of the crystal and the
% prime surface normal in the crystal growth diretion
% refer: Michael A Mastro, Nanotechnology 2021 145205
% https://iopscience.iop.org/article/10.1088/0957-4484/21/14/145205/pdf
function P_piezo=piezo_analytic(theta)

    sin_squared_theta=sin(theta).^2;
    cos_squared_theta=cos(theta).^2;

    sin_2theta=sin(2*theta);
    cos_2theta=cos(2*theta);

    sin_squared_2theta=sin_2theta.^2;
    cos_squared_2theta=cos_2theta.^2;

    sin_fourth_theta=sin(theta).^4;
    cos_fourth_theta=cos(theta).^4;

    a_AlN=3.112;
    a_GaN=3.189;
    a_AlGaN=0.3*a_AlN+0.7*a_GaN;
    a_T=a_GaN;
    a_L=a_AlGaN;

    c_AlN=4.982;
    c_GaN=5.185;
    c_AlGaN=0.3*c_AlN+0.7*c_GaN;
    c_T=c_GaN;
    c_L=c_AlGaN;

    C11_AlN=396;
    C11_GaN=367;
    C11=0.3*C11_AlN+0.7*C11_GaN;

    C12_AlN=137;
    C12_GaN=135;
    C12=0.3*C12_AlN+0.7*C12_GaN;

    C13_AlN=108;
    C13_GaN=103;
    C13=0.3*C13_AlN+0.7*C13_GaN;

    C33_AlN=373;
    C33_GaN=405;
    C33=0.3*C33_AlN+0.7*C33_GaN;

    C44_AlN=116;
    C44_GaN=95;
    C44=0.3*C44_AlN+0.7*C44_GaN;

    e33_AlN=1.55;
    e33_GaN=0.73;
    e33=0.3*e33_AlN+0.7*e33_GaN;

    e31_AlN=-0.58;
    e31_GaN=-0.49;
    e31=0.3*e31_AlN+0.7*e31_GaN;

    e15_AlN=-0.48;
    e15_GaN=-0.40;
    e15=0.3*e15_AlN+0.7*e15_GaN;

    eps_m1=(a_T-a_L)/a_L;

    eps_m2=(a_T*c_T-sqrt((a_L*c_T)^2*cos_squared_theta+(a_T*c_L)^2*sin_squared_theta))./sqrt((a_L*c_T)^2*cos_squared_theta+(a_T*c_L)^2*sin_squared_theta);

    A31=C11*sin_fourth_theta+(0.5*C13+C44)*sin_squared_2theta+C33*cos_fourth_theta;

    A32=(C11*sin_squared_theta+(C13+2*C44)*cos_2theta-C33*cos_squared_theta).*sin_2theta;

    A41=0.5*((C11-C13)*sin_squared_theta+2*C44*cos_2theta+(C13-C33)*cos_squared_theta).*sin_2theta;

    A42=((C11+C33)/2-C13)*sin_squared_2theta+2*C44*cos_squared_2theta;

    B31=C12*sin_squared_theta+C13*cos_squared_theta;

    B32=C13*(sin_fourth_theta+cos_fourth_theta)+((C11+C33)/4-C44)*sin_squared_2theta;

    B41=(C12-C13)/2*sin_2theta;

    B42=0.5*(C11*cos_squared_theta-(C13+2*C44)*cos_2theta-C33*sin_squared_theta).*sin_2theta;

    eps_xx=eps_m1;

    eps_yy=eps_m2;

    eps_zz=((B41*eps_m1+B42.*eps_m2).*A32-(B31*eps_m1+B32.*eps_m2).*A42)./(A31.*A42-A32.*A41);

    eps_yz=((B31*eps_m1+B32.*eps_m2).*A41-(B41*eps_m1+B42.*eps_m2).*A31)./(A31.*A42-A32.*A41);

    P_piezo=eps_xx.*e31*cos(theta)+eps_yy.*(e31*cos(theta).^3+(e33-e15)/2*sin(theta).*sin_2theta)+eps_zz.*((e31+e15)/2*sin(theta).*sin_2theta+e33*cos(theta).^3)+eps_yz.*((e31-e33)*cos(theta).*sin_2theta+e15*sin(theta).*cos_2theta);
end