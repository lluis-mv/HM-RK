function [a, b, c] = GetRungeKutta(method)

substepping = false;
if ( method > 10)
    alfa = floor(method/10);
    method = method-alfa*10;
    substepping = true;
end
    

if ( method == 1)
    % Forward euler
    a = 0;
    b = 1;
    c = 0;
elseif ( method == 2)
    % Heun's method second-order method
    a = [0,0;
        1, 0];
    b = [1/2, 1/2];
    c = [0, 1];
elseif (method == 3)
    % Heun's third-order method
    a = [0, 0, 0;
        1/3, 0, 0;
        0, 2/3, 0];
    b = [1/4, 0, 3/4];
    c = [0, 1/3, 2/3];
elseif (method == 4)
    % Classic fourth-order method
    a = [0,0,0,0;
        1/2,0,0,0;
        0, 1/2, 0, 0;
        0, 0, 1, 0];
    b = 1/6*[1,2,2,1];
    c = [0, 1/2, 1/2, 1];
elseif (method == 5)
    % Fehlberg 5th order method
    a = [0,0,0,0,0,0;
        1/4, 0,0,0,0,0;
        3/32, 9/32, 0,0,0,0;
        1932/2197, -7200/2197, 7296/2197, 0,0,0;
        439/216, -8, 3680/513, -845/4104, 0,0;
        -8/27, 2, -3544/2565, 1859/4104, -11/40,0];
    b = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55];
    c = [0, 1/4, 3/8, 12/13, 1, 1/2];
elseif (method == 6)
    % Dormand and Prince 6th order.
    % https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary/blob/master/Mathematical%20Functions/Differential%20Equations/RungeKStep.m
    a=[0,               0,          0,              0,                  0,                  0,                  0,  0;
        1/10,            0,          0,              0,                  0,                  0,                  0,  0;
        1/36,            5/36,       0,              0,                  0,                  0,                  0,  0;
        10/243,          20/243,     8/81,           0,                  0,                  0,                  0,  0;
        4047/5500,       -18/55,     -4212/1375,     17901/5500,         0,                  0,                  0,  0;
        -5587/4125,      24/55,      9576/1375,      -140049/23375,      38/51,              0,                  0,  0;
        12961/2376,      -35/33,     -160845/5434,   1067565/38896,      -103375/47736,      32875/35568,        0,  0;
        702799/199584,   -1865/2772, -2891375/152152,19332955/1089088,   -5356375/4009824,   2207875/2987712,    0,  0];
    b=[1/12, 0, -216/1235, 6561/12376, 1375/5304, 1375/5928, -5/168, 1/10];
    c=[0; 1/10; 1/6; 2/9; 3/5; 4/5; 1; 1];
elseif (method == 7)
    %Fehlberg's RK7(8)13 formula of [2].
    a=[0,           0,      0,      0,          0,          0,      0,          0,      0,      0,      0,  0,  0;
        2/27,        0,      0,      0,          0,          0,      0,          0,      0,      0,      0,  0,  0;
        1/36,        1/12,   0,      0,          0,          0,      0,          0,      0,      0,      0,  0,  0;
        1/24,        0,      1/8,    0,          0,          0,      0,          0,      0,      0,      0,  0,  0;
        5/12,        0,      -25/16, 25/16,      0,          0,      0,          0,      0,      0,      0,  0,  0;
        1/20,        0,      0,      1/4,        1/5,        0,      0,          0,      0,      0,      0,  0,  0;
        -25/108,     0,      0,      125/108,    -65/27,     125/54, 0,          0,      0,      0,      0,  0,  0;
        31/300,      0,      0,      0,          61/225,     -2/9,   13/900,     0,      0,      0,      0,  0,  0;
        2,           0,      0,      -53/6,      704/45,     -107/9, 67/90,      3,      0,      0,      0,  0,  0;
        -91/108,     0,      0,      23/108,     -976/135    311/54, -19/60,     17/6,   -1/12,  0,      0,  0,  0;
        2383/4100,   0,      0,      -341/164,   4496/1025,  -301/82, 2133/4100, 45/82,  45/164, 18/41,  0,  0,  0;
        3/205,       0,      0,      0,          0,          -6/41,   -3/205,    -3/41,  3/41,   6/41,   0,  0,  0;
        1777/4100,   0,      0,      -341/164,   4496/1025,  -289/82, 2193/4100  51/82,  33/164, 12/41,  0,  1,  0];
    
    b=[41/840, 0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 41/840, 0, 0];
	c=[0; 2/27; 1/9; 1/6; 5/12; 1/2; 5/6; 1/6; 2/3; 1/3; 1; 0; 1];
elseif (method == 8)
    a=[0,                       0,      0,      0,                          0,                      0,                      0,                      0,                      0,                      0,                      0,                      0,  0;
        1/18,                    0,      0,      0,                          0,                      0,                      0,                      0,                      0,                      0,                      0,                      0,  0;
        1/48,                    1/16,   0,      0,                          0,                      0,                      0,                      0,                      0,                      0,                      0,                      0,  0;
        1/32,                    0,      3/32,   0,                          0,                      0,                      0,                      0,                      0,                      0,                      0,                      0,  0;
        5/16,                    0,      -75/64, 75/64,                      0,                      0,                      0,                      0,                      0,                      0,                      0,                      0,  0;
        3/80,                    0,      0,      3/16,                       3/20,                   0,                      0,                      0,                      0,                      0,                      0,                      0,  0;
        29443841/614563906,      0,      0,      77736538/692538347,         -28693883/1125000000,   23124283/1800000000,    0,                      0,                      0,                      0,                      0,                      0,  0;
        16016141/946692911,      0,      0,      61564180/158732637,         22789713/633445777,     545815736/2771057229,   -180193667/1043307555,  0,                      0,                      0,                      0,                      0,  0;
        39632708/573591083,      0,      0,      -433636366/683701615,       -421739975/2616292301,  100302831/723423059,    790204164/839813087,    800635310/3783071287,   0,                      0,                      0,                      0,  0;
        246121993/1340847787,    0,      0,      -37695042795/15268766246,   -309121744/1061227803,  -12992083/490766935,    6005943493/2108947869,  393006217/1396673457,   123872331/1001029789,   0,                      0,                      0,  0;
        -1028468189/846180014,   0,      0,      8478235783/508512852,       1311729495/1432422823,  -10304129995/1701304382,-48777925059/3047939560,15336726248/1032824649, -45442868181/3398467696,3065993473/597172653,   0,                      0,  0;
        185892177/718116043,     0,      0,      -3185094517/667107341,      -477755414/1098053517,  -703635378/230739211,   5731566787/1027545527,  5232866602/850066563,   -4093664535/808688257,  3962137247/1805957418,  65686358/487910083,     0,  0;
        403863854/491063109,     0,      0,      -5068492393/434740067,      -411421997/543043805,   652783627/914296604,    11173962825/925320556,  -13158990841/6184727034,3936647629/1978049680,  -160528059/685178525,   248638103/1413531060,   0,  0];
    b=[14005451/335480064, 0, 0, 0, 0, -59238493/1068277825, 181606767/758867731, 561292985/797845732, -1041891430/1371343529, 760417239/1151165299, 118820643/751138087, -528747749/2220607170, 1/4];
    c=[0; 1/18; 1/12; 1/8; 5/16; 3/8; 59/400; 93/200; 5490023248/9719169821; 13/20; 1201146811/1299019798; 1; 1];
    
end

if (substepping == false)
    return;
end

n = size(a, 1);

a = a/alfa;
b = b/alfa;
A = zeros(alfa*n, alfa*n);
B = zeros(1, alfa*n);

for i = 1:alfa
    ind1 = 1+(i-1)*n:i*n;
    A(ind1, ind1) = a;
    B(ind1) = b;
    
    for j = 1:i-1
        ind2 = 1 +(j-1)*n:j*n;
        for k = ind1
            A(k,ind2 ) = b;
        end
    end
    
end


a = A;
b = B;

