%% Performance of Slotted ALOHA on infrared WBANs
% This simulation tries to describe the performance of Slotted ALOHA
% when applied as a MAC protocol for a WBAN system that consists of
% 6 nodes. For the communication medium it is assumed that we use the
% infrared part of the electromagnetic spectrum.

% Authors:
% Giachoudis Christos, christos.giachoudis@fresnel.fr
% Vasilis Papanikolaou,
% Konstantinos Rallis,


%% Equations

% Channel Statistical Model (Gamma distribution)
% To describe the channel DC gain we use the gamma distribution
x = 1:10;
a = 13.79;
b = 0.04;
f_x = gampdf(x, a, b);

% Average throughput of the network Rk_hat
% k is the number of the kth node
% Rk_power is the average rate of the node
% q(k) is the probability of channel access of the kth node
function Rk_hat = avThrouput (k, Rk_power, q)
    temp1 = Rk_power * q(k);
    temp2 = 1;
    for i = 1:k
        if i ~= k
            temp2 = temp2 * (1 - q(i));
        end
    end
    Rk_hat = temp1 * temp2;
end


% Average rate of the kth node Rk_power

function Rk_power = avRate (Rk, k, Pk)
    temp1 = ;
    temp2 = 1 - temp1;
    Rk_power = Rk * temp2;
end
