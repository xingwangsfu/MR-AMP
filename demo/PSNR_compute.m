PSNR_HR = 20*log10(255/sqrt(mean((alpha_rec(:)-alpha_true(:)).^2)))
PSNR_L2H = 20*log10(255/sqrt(mean((alpha_L2H(:)-alpha_true(:)).^2)))
alpha_extend = zeros(128,128);
alpha_LR_mtx = reshape(alpha_rec_LR,64,64);
alpha_extend(1:64,1:64) = alpha_LR_mtx;
PSNR_L2H_extend = 20*log10(255/sqrt(mean((alpha_extend(:)-alpha_true(:)).^2)))