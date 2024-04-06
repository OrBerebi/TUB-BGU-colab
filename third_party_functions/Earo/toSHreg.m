        function obj = toSHreg(obj, orderN, dim,regFlag,regParam,cond_flag)
            % function obj = toSHreg(obj, orderN, dim,regFlag)
            %
            % transform data to the spherical harmonics (SH) domain from
            % the space domain using regularization
            %
            %  orderN    - transformation order
            %  dim       - transformation dimension ('MIC' or 'SRC')
            %  regFlag   - regularization flag: 1 SVD; 2- DL; 3- robust; 0- conventional SFT
            %  regParam  - regularization parameter
            if nargin<4,
                regFlag=0;                
            end
            if nargin<6,
                cond_flag=0;
            end
            
            %%
            if strcmp(obj.dataDomain{2}, 'SPACE')
                if strcmp(dim,'SRC')
                    
                    if ~obj.shutUp, fprintf('Transforming to SH over source grid, quadrature = %s\n',obj.sourceGrid.quadType); end;
                    el = obj.sourceGrid.elevation;
                    az = obj.sourceGrid.azimuth;
                    Yq=(sh2(orderN, el, az)).';
                    switch regFlag
                        case 1
                            if nargin<5,
                                regParam=9;
                            end
                            Nl=regParam;
                            [U,S,V]=svd(Yq);   %Thus, Yq=U*S*V' and pinv(Yq)=V*inv(S'*S)*S'*U'~V*inv(S'*S)*S'*U' 
                            S2=diag(abs(S).^2);
                            invS2=zeros((orderN+1)^2);
                            invS2(1:(Nl+1)^2,1:(Nl+1)^2)=diag(S2(1:(Nl+1)^2).^-1);
                        Yp=V*invS2*S'*U'; %By zeroing the low SV approximate the expression: pinv(Yq)=V*inv(S'*S)*S'*U'    
                        case 2
                            if nargin<5,
                                regParam=0.001;
                            end
                            
                            [U,S,V]=svd(Yq);   %Thus, Yq=U*S*V' and pinv(Yq)=V*inv(S'*S)*S'*U'~V*inv(S'*S)*S'*U'
                            dl=(regParam*S(1,1))^2; % The diagonal loading is set according to the highest desired SV (10% of its value)
                            S2=(S'*S)+dl*eye((orderN+1)^2);
                            Yp=V*inv(S2)*S'*U'; %By zeroing the low SV approximate the expression: pinv(Yq)=V*inv(S'*S)*S'*U'
                            if cond_flag,
                                display(['Cond of matrix Y: cond(Y)=',num2str(round(cond(Yq))),', and after diagonal loading cond(Yd)=',num2str(round(cond(Yp)))]);
                            end
                        case 0
                            Yp=pinv(Yq);
                            if cond_flag,
                                display(['Cond of matrix Y: cond(Y)=',num2str(round(cond(Yq))),', No diagonal loading was performed']);
                            end
                    end
                    
                    for nn=1:size(obj.data,3)       % iterate through mics
                        tmp(:,:,nn) = Yp*squeeze(obj.data(:,:,nn));
                    end
                    obj.data = tmp;
                    obj.orderN = orderN;
                    obj.dataDomain{2} = 'SH';
                    
                    
                elseif strcmp(dim,'MIC')
                    % Perhaps change this to transformation using B?
                    if ~obj.shutUp
                        fprintf('Transforming to SH over mic grid, quadrature = %s\n',obj.micGrid.quadType);
                        if obj.isRobot
                          fprintf('-- This transformation assumes a spherical distribution of mics, but isRobot is true! --\n');
                        end
                    end
                    el = obj.micGrid.elevation;
                    az = obj.micGrid.azimuth;
                    Yp=pinv(sh2(orderN, el, az));    
                    for nn=1:size(obj.data,1)       % iterate through sources
                        tmp(nn,:,:) = (squeeze(obj.data(nn,:,:))*Yp);
                    end
                    obj.data = tmp;
                    obj.orderN = orderN;
                    obj.dataDomain{2} = 'SH';
                    
                else
                    if ~obj.shutUp, disp ('Invalid transformation dimension. Please specify "MIC" or "SRC"'); end;
                end
            else
                if ~obj.shutUp, disp('Warning! Data already in SH domain'); end;
            end
        end

