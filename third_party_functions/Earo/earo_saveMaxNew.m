function earo_saveMaxNew(eobj, filename, hpstruct)

    % saveSSR(eobj, filename, norml, [hpobj])
    %
    % based on earo_saveSSR by Jonathan Sheaafer
    % edited by Zamir Ben-Hur 1.3.2016
    %
    
  %  filenameLeft = [filename,'Left.wav'];
   % filenameRight = [filename,'Right.wav'];
   % nfft = 2^12;
    % resapmle to 44100 Hz;
   % newFs = 44100;
   % [P,Q] = rat(newFs/eobj.fs);
    
    if (strcmp(eobj.type,'BRIR')  || strcmp(eobj.type,'HRIR')) && eobj.nData == 360

        if ~eobj.shutUp && nargin>=4, fprintf('Applying headphone compensation, %s\n',hpstruct.hpName); end;

        irArrayLeft = zeros(size(eobj.data,2), 360);
        irArrayRight = zeros(size(eobj.data,2), 360);
        inc = 1;

        for i = 1:360
            if ~mod(i,10)
                fprintf('|');
            end
            
            IR = squeeze(eobj.data(i,:,:));

            if nargin>=3        % Headphone compensation
                if hpstruct.fs == eobj.fs
                    irArrayLeft(:,inc) = miro_fastConv(IR(:,1),hpstruct.minPhase);
                  
                    %Right ear
                    irArrayRight(:,inc) = miro_fastConv(IR(:,2),hpstruct.minPhase);
                    
                    inc = inc+1;       
                else
                    fprintf ('Headphone filter and earo object sample rate mismatch!\nDoing nothing.\n')
                end
            else
                irArrayLeft(:,inc) = IR(:,1);
                    
                %Right ear
                irArrayRight(:,inc) = IR(:,2);
                    
                inc = inc+1;
            end

        end

        fprintf('\n\n');
        

        irArrayLeft = 0.99*irArrayLeft/max(max(abs(irArrayLeft)));
        irArrayRight = 0.99*irArrayRight/max(max(abs(irArrayRight)));
        
        dataforMaxLeft = zeros(size(irArrayLeft));
        dataforMaxRight = zeros(size(irArrayRight));
        
        for i = 1:181
            dataforMaxLeft(:,182-i) = irArrayLeft(:,i);
            dataforMaxRight(:,182-i) = irArrayRight(:,i);
        end
        for i = 2:180
            dataforMaxLeft(:,362-i) = irArrayLeft(:,i+180);
            dataforMaxRight(:,362-i) = irArrayRight(:,i+180);
        end
        
        mkdir(filename);
        for i=1:360
            wavwrite(dataforMaxLeft(1:end,i), eobj.fs, 16, [filename,'Left',num2str(i-1),'.wav'])
            wavwrite(dataforMaxRight(1:end,i),eobj.fs, 16, [filename,'Right',num2str(i-1),'.wav'])
        end
        disp(['Max file sucessfully generated: ', filename])
    else
        if ~obj.shutUp
            fprintf('Failed: Max export requires circular HRIR or BRIR objects.\n');
        end
        return
    end

end

function ab = miro_fastConv(a,b)

% Internal use only

NFFT = size(a,1);%+size(b,1)-1;
A    = fft(a,NFFT);
B    = fft(b,NFFT);
AB   = A.*B;
ab   = ifft(AB);

end