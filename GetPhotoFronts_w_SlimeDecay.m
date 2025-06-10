function [ZeroContouryCell,TimeVec,Xvec,Yvec,Final_C, Final_S, OpStruct, OpFname] = GetPhotoFronts_w_SlimeDecay(Cmean_nondim,fc_nondim, tau_nondim, graphout)

    %Ritwika VPS; Jun 2025 (Modified from Tristan Ursell's code (see below))
    % This function gets bacterial fronts from 2d phototaxis simulations with slime decay (because that is the equivalent model to teh 1d model). Note that this realisation of the 
    % 2d model has no reproduction and as such, the total number of bacterial cells remain constant.
    
    %-------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Tristan Ursell; December 2010
    % Reaction-Diffusion for Cyanobacteria motility where cell concentration and slime concentration is in volume per area
    %-------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    % Unless otherwise stated, the source for the simulation framework and parameter values is the paper titled: Motility Enhancement through Surface Modification Is Sufficient 
    % for Cyanobacterial Community Organization during Phototaxis, Ursell et al, 2013, PLOS Comp Bio (referred to as Motility enhancement paper in the rest of this script).
    
    % Inputs: - Cmean_nondim: non-dimensionalised value of the initial mean concentration that we start the simulations with
    %                         (denoted by C_tot_tilde in the Motility enhancement paper, and expressed as a power of 2, to easily match with Figs 6 and 9 in the Motility 
    %                         enhancement paper; see below):
    %                               eg. Cmean_nondim = 2^(-1) for simulation profile type 3 in the phase diagram in Fig. 9; units prior to nondimensionalisation is um. 
    %                               See TU's original script(s) for preliminary estimates of initial concentration
    %         - fc_nondim: non-dimensionalised value of the bias forced (denoted by beta_tilde in the paper, with the pre-nondiensionalisation vector version with units of 1/um
    %                      denoted by b in vector notation, and the magnitude of b denoted as beta; expressed as a power of 2, to easily match with Figs 6 and 9 in the Motility 
    %                      enhancement paper)
    %                               eg. fc_nondim = 4 (for Fig. 9, mentioned in pg 10 of the paper).
    %         - tau_nondim: non-dimensionalised slime decay constant (represented by tau prior to non-dimensionalisation and by tau_tilde after non-dimensionalisation in the paper;
    %                       expressed as a power of 2, to easily match with Fig 9 in the Motility enhancement paper)
    %                               eg. tau_nondim = 2^(0), for simulation profile type 3 in the phase diagram in Fig. 9; units os seconds prior to nondimensionalisation
    % 
    % Outputs: ZeroContouryCell: The cell array with the fronts (estimated as the Y point at which the concentration falls to zero for that Y transect) at each time step. In um, 
    %                            because Xvec and Yvec are in um.
    %          TimeVec: Time points (in s; converted from simulation step index) corresponding to each front in ZeroContouryCell
    %          Xvec, Yvec: The vector of X and Y grid points (um)
    %          Final_C: Final conc. profile (at final time point achieved in simulation, whether simultation ends before the full simulation time (600000 s) because fronts reach 
    %                   the end of simulation area or if the simulation proceeds till the end of simulation time.
    %          Final_S: Similarly as above, final slime profile
    %          OpStruct: output structure with various simulation parameters and outputs stored
    %          OpFname: file name for the output structure
    %          graphout: graph out logical (true outputs a graph).
    
    
    fout = 0; %file out logical
    circ = 0; %whether intial concentration is initialised as circular or not, logical 
    WdthToHtRatio_Lt = 1; %set the ratio of the width of a finger (feature) to the height. If there is a feature in the front that has a lower ratio than this,
    %the front is assessed to be non-linear. 
    
    %-------------------------------------------------------------------------------------------------------------------------------------------------------------
    %Define fixed parameters
    %-------------------------------------------------------------------------------------------------------------------------------------------------------------
    k_s = 0.0005; %slime production rate (1/s) (converted from ~ 0.03 per min; source: Motility enhancement paper, see Methods)
    m0 = 0.25; %saturated mobility (um^2/s) (source: Motility enhancement paper, see Methods)
    
    T_nat = 1/k_s; %natural time scale (s) (source: Motility enhancement paper, see Model description under Results)
    X_nat = sqrt(m0/k_s); %natural length scale (um) (source: Motility enhancement paper, see Model description under Results)
    
    sigma = 0.25; %Slime saturation depth, um.
    % I cannot find a source for this number in the paper, but the simulations seem to work when the value of sigma is set to 0.25. Also, note that in Tristan's code, this value 
    % has been rescaled by the natural length scale (see commented out code below). However, while this seems to reproduce simulation results for at least Cmean_nondim = 2^-1, 
    % scaled fc = 4, and scaled slime decay rate = 2^1, the mobility field is very high throughout the simulation area at the end of the simulation because the saturation depth 
    % is very low (~ 0.01).
    % sigma=sig/X_nat;
    
    %Universal time and distance steps, for the simulation.
    dt = 0.01; %Universal time step (dimensionless). This is how many units of the natural time scale makes up one simulation time step.
    dx = 0.25; %universal distance step (dimensionless). This is how many units of the natural length scale makes up one simulation grid square size.
    
    %Get the number of time steps
    tact = 600000; %Actual time simulation, in seconds (source: Motility enhancement paper, see Methods)
    N = round(tact/T_nat/dt); %round to nearest 10
    
    %Set up phototactic bias
    theta = (pi/2)*ones(1,N); %light direction
    F0 = fc_nondim*[cos(theta);sin(theta)];
    
    %Number of x and y grid points
    lx = 2000; ly = 2*lx; %um (2 x 4 mm; source: Motility enhancement paper, see Methods)
    nx = 10*ceil(lx/X_nat/dx/10); ny = 10*ceil(ly/X_nat/dx/10);  %round to nearest 10
    
    
    %-------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Set up display items and related logistics
    %-------------------------------------------------------------------------------------------------------------------------------------------------------------
    %small parameter (to check for simulation stability)
    dsmall=dt/dx^2;
    
    %Frequency of outputting concentration profiles as a heat map
    savestep=round(N/500);
    
    %Set up output structure
    OpStruct.T_nat = T_nat; OpStruct.X_nat = X_nat;
    OpStruct.sigma = sigma;
    OpStruct.fc_nondim = fc_nondim; OpStruct.Cmean_nondim = Cmean_nondim; OpStruct.tau_nondim = tau_nondim; 
    OpStruct.dx = dx; OpStruct.dt = dx;
    
    %Set up base file name
    tdata = datetime; %get date and time info to save files if needed
    fdate = [date '_' num2str(tdata.Hour) '-' num2str(tdata.Minute) '-' num2str(tdata.Second)];
    OpFname = ['2dModelSims_w_SlimeDecay_' fdate '.mat'];
    
    disp(' ')
    disp(['The non-dimensionalised mean initial cell concentration, Cmean_nondim = ' num2str(Cmean_nondim) '.'])
    disp(['The force instability length scale is ' num2str(X_nat/fc_nondim) ' um, fc_nondim = ' num2str(fc_nondim) '.'])
    disp(['Frame frequency: ' num2str(savestep) ' steps, ' num2str(dt*savestep) ' T_nat, ' num2str(dt*savestep*T_nat) ' secs.'])
    %disp(['The length scale is ' num2str(X_nat) ' um.'])
    %disp(['The time scale is ' num2str(T_nat) ' s.'])
    %disp(['The total simulation time is ' num2str(N*T_nat*dt) ' s, with N = ' num2str(N) '.'])
    %disp(['The box size is ' num2str(dx*nx*X_nat) ' x ' num2str(dx*ny*X_nat) ' um, grid size ' num2str(nx) ' x ' num2str(ny) '.'])
    %disp(['The slime depth constant is ' num2str(sigma*X_nat) ' um, sigma = ' num2str(sigma) '.'])
    
    if dsmall > 0.25
        disp(['This computation might be unstable, the small parameter is ' num2str(dsmall) '.'])
        disp(' ')
        input(' ','s')
    else
        disp(['This computation should be stable, the small parameter is ' num2str(dsmall) '.'])
        disp(' ')
        %input(' ','s')
    end
    
    
    %-------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Set up for computation
    %-------------------------------------------------------------------------------------------------------------------------------------------------------------
    %Initialize matrices
    C = zeros(ny,nx); %Bacterial concentration field (spanning the entire simulation grid), %um
    S = 0/X_nat*ones(ny,nx); %Slime concentration field (spanning the entire simulation grid), %um
    M = zeros(ny,nx); %Mobility field (spanning the entire simulation grid)
    dcdx = zeros(ny,nx); %dC/dx
    dcdy = zeros(ny,nx); %dC/dy
    
    %Define position vectors/matrices
    Xvec = X_nat*(dx/2:dx:(nx-1)*dx+dx/2);
    Yvec = X_nat*(dx/2:dx:(ny-1)*dx+dx/2);
    Xpos = zeros(ny,nx);
    Ypos = zeros(ny,nx);
    for i = 1:nx
        Ypos(:,i) = Yvec;
    end
    for i = 1:ny
        Xpos(i,:) = Xvec;
    end
    
    %Setup initial conditions
    thick = round(ny/5);
    C(1:thick,:) = 2*Cmean_nondim*rand(thick,nx); %top rectangle
    
    %Top circle (if circular initial condition)
    if circ == 1
        R = nx/2.1;
        Rx = nx/2; Ry = 1.03*R;
        for i = 1:ny
            for j = 1:nx
                Dmask(i,j) = sqrt((i-Ry)^2+(j-Rx)^2);
            end
        end
        Cmask = Dmask <= R;
        Cini = 2*Cmean_nondim*rand(ny,nx);
        C = Cini.*Cmask;
    end
    
    %Setting up heat map colorbar limits
    C0 = sum(sum(C));
    Clow = 0; Chigh = 1.25*C0/nx/2;
    Slow = 0; Shigh = sigma*5;
    Shigh = 0.14; Chigh = 0.065;
    
    %Open figure if figure display logical is true
    if graphout == 1
        h1 = figure;
        colormap(hot)
    else
        h = waitbar(0,'Calculate cell motility...');
    end
    
    
    %-------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Do the simulation
    %-------------------------------------------------------------------------------------------------------------------------------------------------------------                                                                  
    tic
    
    Ctr = 0;
    for i = 1:N %for N time steps; note that EVERYTHING is non-dimensionalised
        
        S = S*(1 - (dt/tau_nondim)) + C*dt; %calculate the slime field from the concentration field
        %C=C*(1+1/tdoub); %perform growth
        
        M = 1 - exp(-S); %calculate the mobility field from the slime field
        %M=1-exp(-S/sigma); %Previously, there was this extra non-dimensionalisation in Tristan's code, which is not necessary, since everything is already non-dimensionalised 
                                                                        
        Fx = F0(1,i)*ones(ny,nx-1); Fy = F0(2,i)*ones(ny-1,nx); %Specify force field
        
        %Explicitly calculate directional first derivatives, non-symmetric
        dcdx = (C(:,2:end) - C(:,1:end-1))/dx; 
        dcdy = (C(2:end,:) - C(1:end-1,:))/dx;
        
        %Calculate neighbor mean mobility field (with one less element in each direction as 'C')
        Mx = 1/2*(M(:,1:end-1) + M(:,2:end));
        My = 1/2*(M(1:end-1,:) + M(2:end,:));
        
        %Calculate neighbor mean concentration field (with one less element in each direction as 'C')
        Cx = 1/2*(C(:,1:end-1) + C(:,2:end));
        Cy = 1/2*(C(1:end-1,:) + C(2:end,:));
        
        %Calculate flux components (with one less element in each direction as 'C')
        Jx = Mx.*(-Cx.*Fx + dcdx);
        Jy = My.*(-Cy.*Fy + dcdy);
        
        %Update concentration field: x component
        dCx = zeros(ny,nx); 
        dCx(:,1) = Jx(:,1);
        dCx(:,end) = -Jx(:,end);
        dCx(:,2:end-1) = Jx(:,2:end) - Jx(:,1:end-1);
    
        %Update concentration field: y component
        dCy = zeros(ny,nx); 
        dCy(1,:) = Jy(1,:);
        dCy(end,:) = -Jy(end,:);
        dCy(2:end-1,:) = Jy(2:end,:) - Jy(1:end-1,:);
        
        %Update concentration field;
        C = C + dt/dx*(dCx + dCy);
        
        %End detection if bacterial cells have reached the end (Y direction) of simulation grid.
        if sum(C(end,:))>0
            break
        end
        
        %---------------------------------------------------------------------
        %This block was previously in the savestep loop, but moving here for more data, at more time points, to do FFT.
        Ctr = Ctr + 1; %Update Counter
            
        %Store concentration profiles
        if ~isempty(C) %continue only if non-empty
            for j = 1:size(C,2) %each y transect is a column, so we loop over number of columns
                % Previous attempts to find boundaries based on max negative slope or even the last point where there is a non-zero slope. Currently, I am choosing to simply
                % identify the cell boundary as the last Y point where the cell concentration drops to zero. This introduces some 'smearing' in that the boundary is not ultra-precise
                % because the fronts are not such that there is a sudden drop in concetration to 0. However, this might work well-enough for detecting wave numbers. If it does not
                % work well, we will have to detect with better precision, and the gradient might be a good measure for that (as shown below in the commented out portion).
                % Cdiff(:,j) = diff(C(:,j)); %get gradient (since the denominator is the same, just take difference)
                % %[mm,in] = max(-Cdiff(:,j));
                % NonZeroInds = find(Cdiff(:,j));
                % contourymax(j) = Yvec(max(NonZeroInds)+1); %find contour based on max slope
    
                NonZeroInds = find(C(:,j)); %find indices of non-zero concentration
                if max(NonZeroInds) >= size(C,1) %if this is the max Y axis value, then the boundary is at the Y edge of the simulation
                    ZeroContourY(j) = Yvec(max(NonZeroInds));
                else %if not
                    ZeroContourY(j) = Yvec(max(NonZeroInds)+1); %pick the next Y value as the boundary for that Y transect
                end
                ZeroContourY(j) = Yvec(max(NonZeroInds)+1); 
            end
            %CdiffCell{counter} = Cdiff; %contourymaxCell{Ctr} = contourymax;
            ZeroContouryCell{Ctr} = ZeroContourY;
            NonLinOrNo(Ctr) = FrontNonLinTest_PkWdthsAndHts(Xvec, ZeroContourY, WdthToHtRatio_Lt);
            clear Cdiff ZeroContourY %contourymax
        end
    
        TimeVec(Ctr) = i*T_nat*dt; %Get the (unscaled) time vector for fitting for growth rates after Fourier analysis
        
        %---------------------------------------------------------------------
    
        %Save/produce plots every 'savestep'
        if mod(i,savestep) == 1
    
            %Plots
            if graphout == 1
                subplot(1,2,1, 'Parent', h1)
                imagesc(Xvec,Yvec,C); %imagesc(Xvec,Yvec,C,[Clow,Chigh])
    
                hold on
                dn = 10;
                quiver(Xpos(1:dn:end-1,1:dn:end-1),Ypos(1:dn:end-1,1:dn:end-1),-Jx(1:dn:end-1,1:dn:end),-Jy(1:dn:end,1:dn:end-1))
                hold off
    
                xlabel('x (um)'); ylabel('y (um)')
                title(['Concentration, \Delta =  ' num2str(sum(sum(C))/C0) ' , t = ' num2str(i*dt) '/' num2str(N*dt) ', T_nat (s) = ' num2str(T_nat)])
                axis equal; axis tight; box on
                
                subplot(1,2,2, 'Parent', h1)
                imagesc(Xvec,Yvec,S) %imagesc(Xvec,Yvec,S,[Slow,Shigh])
                xlabel('x (um)'); ylabel('y (um)')
                title(['Slime, fc_nondim = ' num2str(fc_nondim) ', \sigma = ' num2str(sigma) ', <C> = ' num2str(Cmean_nondim) ', X_nat (um) = ' num2str(X_nat)])
                axis equal; axis tight; box on
                colormap(hot)
                pause(0.1)
            else
                waitbar(i/N,h)
            end
            
            %Saving file
            if fout == 1 %write output image
                imwrite(mat2gray(C,[Clow,Chigh]),[fname '_C.tiff'],'Compression','none','WriteMode','append');
                imwrite(mat2gray(S,[Slow,Shigh]),[fname '_S.tiff'],'Compression','none','WriteMode','append');
            end
        end
    end
    
    %Get the final concentration and slime profiles
    Final_C = C; Final_S = S;

    if graphout == 0
        close(h)
    elseif graphout == 1
        close(h1) %close figure
    end
    
    %Get computation time
    tel = toc;
    disp(['Computation took: ' num2str(tel) ' s'])
    
    %Set up the rest of the outputs to save
    OpStruct.FingerContour = ZeroContouryCell;
    OpStruct.Time = TimeVec;
    OpStruct.Xvec = Xvec;
    OpStruct.Yvec = Yvec;
    
    %Waiting for user input before exiting
    %input('Please check to make sure that finger formation occured for this combination of parameters',"s")

end


