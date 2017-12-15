function mask  = load_mask( mask_type, pars )
%
%   load_mask: 
%
%       mask_type: 'random', 'radial', 'slice'
%

m = pars.image_size(1); n = pars.image_size(2);
% Generate mask by type and parameters
switch mask_type 
    case 'random'
        OMEGA = RandMask_rect(double(m/pars.d),double(n/pars.d),m,n);
        mask = RandMask_InverseTransfer(OMEGA,m,n);
    case 'radial'
        OMEGA = fftshift(MRImask(m, pars.line_num));
        mask = RandMask_InverseTransfer(OMEGA,m,n);
    case 'cartesian'
        pdf = genPDF([1,m], 6, 1/pars.d, 2, 0.1, 0);
        mask = fftshift(fftshift(repmat(genSampling_symmetric(pdf,10,1),[m,1])));
        if abs(m-n)>0
            mask = imresize(mask,[m,n]);
        end
    otherwise
        error('Unsupported mask type');
end

% Construct Image mask by fftshift

% OMEGA = find(fftshift(mask)==1); % Currently no use

% Using central window mask for e.g. GRAPPA
if isfield(pars, 'central_window') 
    cw_x = pars.central_window(1);
    cw_y = pars.central_window(2);
    mask(m/2-cw_x/2+1:m/2+cw_x/2,n/2-cw_y/2+1:n/2+cw_y/2) = ones(cw_x, cw_y);
end

end

