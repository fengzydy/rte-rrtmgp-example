
function [olr]=execute_RRTMGP_lw_input(lon,lat,p_lev,p_lay,sfc_t,t_lay,q_lay,o3);
% output: the name of output nc file
% Vertical layers are inputted from TOA to ground
% p_lev: level-by-level pressure [Pa]
% p_lay: layer-by-layer pressure [Pa]
% sfc_t: surface temperature [K]
% t_lay: layer-by-layer temperatures [K]
% q_lay: layer-by-layer water vapor mixing ratio [kg/kg]
% o3: ozone mass mixing ratio [kg/kg]
    write_input_atmos('test.nc',lon,lat,p_lev,p_lay,sfc_t,t_lay,q_lay,o3);
    system('./run_clear_lw.sh');
    tmp = ncread('test.nc','lw_flux_up');
    olr = tmp(:,1);
    nlon=length(lon);
    nlat=length(lat);
    olr           = reshape(olr,[nlon,nlat]);
end