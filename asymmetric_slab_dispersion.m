function asymmetric_slab_dispersion()
% Reconstructed from the supplied slab-waveguide report.
% Numerical dispersion analysis of TE and TM modes.

clc;
close all;

disp('ASYMMETRIC SLAB WAVEGUIDE DISPERSION');

% User inputs
n1 = input('Enter core refractive index n1: ');
n_top = input('Enter TOP cladding index n_top: ');
n_bottom = input('Enter BOTTOM cladding index n_bottom: ');
d = input('Enter slab thickness d (meters): ');

params = [n1, n_top, n_bottom, d];
if any(~isfinite(params)) || any(params <= 0)
    error('All inputs must be positive finite values.');
end

if ~(n1 > n_top && n1 > n_bottom)
    disp('No guided modes: core index must exceed both claddings.');
    return;
end

m_max = 4;
V = linspace(0.01, 10, 250);
n_ref = min(n_top, n_bottom);

b_TE = NaN(m_max+1, numel(V));
b_TM = NaN(m_max+1, numel(V));

for k = 1:numel(V)
    Vnow = V(k);
    k0 = Vnow / (d * sqrt(n1^2 - n_ref^2));

    % Bound modes must lie above both cladding light lines.
    beta_min = k0 * max(n_top, n_bottom) * (1 + 1e-8);
    beta_max = k0 * n1 * (1 - 1e-8);

    if beta_min >= beta_max
        continue;
    end

    beta_scan = linspace(beta_min, beta_max, 4000);
    Fte = arrayfun(@(b) TE_disp(b,k0,d,n1,n_top,n_bottom), beta_scan);
    Ftm = arrayfun(@(b) TM_disp(b,k0,d,n1,n_top,n_bottom), beta_scan);

    for m = 0:m_max
        % TE roots
        G = Fte - m*pi;
        idx = find(isfinite(G(1:end-1)) & isfinite(G(2:end)) & ...
                   G(1:end-1).*G(2:end) < 0);

        if ~isempty(idx)
            j = idx(1);
            try
                beta_root = fzero(@(b) ...
                    TE_disp(b,k0,d,n1,n_top,n_bottom)-m*pi, ...
                    [beta_scan(j), beta_scan(j+1)]);
                neff = beta_root/k0;
                b_TE(m+1,k) = ...
                    (neff^2-n_ref^2)/(n1^2-n_ref^2);
            catch
            end
        end

        % TM roots
        G = Ftm - m*pi;
        idx = find(isfinite(G(1:end-1)) & isfinite(G(2:end)) & ...
                   G(1:end-1).*G(2:end) < 0);

        if ~isempty(idx)
            j = idx(1);
            try
                beta_root = fzero(@(b) ...
                    TM_disp(b,k0,d,n1,n_top,n_bottom)-m*pi, ...
                    [beta_scan(j), beta_scan(j+1)]);
                neff = beta_root/k0;
                b_TM(m+1,k) = ...
                    (neff^2-n_ref^2)/(n1^2-n_ref^2);
            catch
            end
        end
    end
end

% Plot normalized propagation constant
figure;
hold on;
grid on;
box on;

for m = 0:m_max
    if any(isfinite(b_TE(m+1,:)))
        plot(V,b_TE(m+1,:),'b-','LineWidth',1.2);
    end
    if any(isfinite(b_TM(m+1,:)))
        plot(V,b_TM(m+1,:),'r--','LineWidth',1.2);
    end
end

title('Normalized Propagation Constant b vs V');
xlabel('Normalized Frequency V');
ylabel('Normalized Propagation Constant b');
ylim([0 1.05]);
legend('TE modes','TM modes');
hold off;

end

function F = TE_disp(beta,k0,d,n1,n_top,n_bottom)
h2 = n1^2*k0^2-beta^2;
a2 = beta^2-n_top^2*k0^2;
g2 = beta^2-n_bottom^2*k0^2;

if h2 <= 0 || a2 <= 0 || g2 <= 0
    F = NaN;
    return;
end

h = sqrt(h2);
alpha = sqrt(a2);
gamma = sqrt(g2);

F = h*d-(atan(alpha/h)+atan(gamma/h));
end

function F = TM_disp(beta,k0,d,n1,n_top,n_bottom)
h2 = n1^2*k0^2-beta^2;
a2 = beta^2-n_top^2*k0^2;
g2 = beta^2-n_bottom^2*k0^2;

if h2 <= 0 || a2 <= 0 || g2 <= 0
    F = NaN;
    return;
end

h = sqrt(h2);
alpha = sqrt(a2);
gamma = sqrt(g2);

F = h*d-( ...
    atan((n1^2/n_top^2)*(alpha/h)) + ...
    atan((n1^2/n_bottom^2)*(gamma/h)));
end
