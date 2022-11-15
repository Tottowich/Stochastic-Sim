clear all
close all

f = @(x) exp(-(1/2).*x.^2).*(1./sqrt(2*pi));
g = @(x) (1/pi)*1./(1+x.^2);

c = sqrt(2*pi./(exp(1)));

N = 100000;
X = zeros(1,N);
% plot densities f(x), g(x) and function c*g(x)
figure(1)
t=-20:0.01:20;
plot(t,f(t),t,g(t),t,c*g(t))
legend("f","g","c*g")
grid on
for i = 1:N
    while true
        u = rand;
        y = tan(pi*(u-1/2));
        u = rand;
        if c*g(y)*u<f(y) % Acceptance rejection method.
            x=y;
            break
        end
        
    end
    X(i) = x;
end
figure(2)
histogram(X,'Normalization','pdf')
hold on
plot(t,f(t),'LineWidth',3)
grid on
hold off