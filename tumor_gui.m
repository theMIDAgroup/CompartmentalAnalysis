function tumor_gui

global analysis_gui
%-------------------------------------------------------------------------%
% whitecolor=[1,1,1];
startstop=false;
critval=1e-4;
relerr=Inf;
nitmax=500;
%V=; %Vb blood fraction
[glnodes,glweights]=gauss_legendre(4);
alpha=[1,1];
%-------------------------------------------------------------------------%
t0=0;
%-------------------------------------------------------------------------%
[t,Ca,Cdata]=deal(0);
[Kfbx,Kmfx,Kbfx,Kfmx,Cx]=deal(0);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
h = figure('Visible','off','Position',[200,200,690,450],...
    'MenuBar','none','Toolbar','figure');
%-------------------------------------------------------------------------%

% TUMOR
%-------------------------------------------------------------------------%
htumorpanel = uipanel('Title','Tumor',...
    'TitlePosition','centertop',...
    'Position',[460/690,63/450,220/690,360/450]);
%-------------------------------------------------------------------------%
htumorig = uicontrol('Style','text','String','Initial guess',...
    'Position',[5,315,110,20],...
    'HorizontalAlignment','Left','Parent',htumorpanel);
%-------------------------------------------------------------------------%
htumorigrand = uicontrol('Style','pushbutton','String','Random',...
    'Position',[110,315,105,25],...
    'HorizontalAlignment','Right','Parent',htumorpanel,...
    'Callback',{@tumorigrand_Callback});
%-------------------------------------------------------------------------%
hkmftxt = uicontrol('Style','text','String','Kmf=',...
    'Position',[0,280,35,20],...
    'Parent',htumorpanel);
hkmf = uicontrol('Style','edit','String','',...
    'Position',[35,280,70,25],...
    'Parent',htumorpanel);
%-------------------------------------------------------------------------%
hkfmtxt = uicontrol('Style','text','String','Kfm=',...
    'Position',[0,250,35,20],...
    'Parent',htumorpanel);
hkfm = uicontrol('Style','edit','String','',...
    'Position',[35,250,70,25],...
    'Parent',htumorpanel);
%-------------------------------------------------------------------------%
hkbftxt = uicontrol('Style','text','String','Kbf=',...
    'Position',[110,280,35,20],...
    'Parent',htumorpanel);
hkbf = uicontrol('Style','edit','String','',...
    'Position',[145,280,70,25],...
    'Parent',htumorpanel);
%-------------------------------------------------------------------------%
hkfbtxt = uicontrol('Style','text','String','Kfb=',...
    'Position',[110,250,35,20],...
    'Parent',htumorpanel);
hkfb = uicontrol('Style','edit','String','',...
    'Position',[145,250,70,25],...
    'Parent',htumorpanel);
%-------------------------------------------------------------------------%
hregtxttit = uicontrol('Style','text','String','Regularization parameter',...
    'Position',[5,215,210,20],...
    'Parent',htumorpanel);
hregtxt = uicontrol('Style','text','String','r=',...
    'Position',[70,190,30,20],...
    'Parent',htumorpanel);
hreg = uicontrol('Style','edit','String','1e4',...
    'Position',[100,190,70,25],...
    'Parent',htumorpanel);
%-------------------------------------------------------------------------%
htumorstartstop = uicontrol('Style','pushbutton','String',...
    'Start/Stop iterations','Position',[10,155,200,25],...
    'Parent',htumorpanel,...
    'Callback',{@tumorstartstop_Callback});
%-------------------------------------------------------------------------%
htumornittxt = uicontrol('Style','text','String','Iteration',...
    'Position',[15,130,105,20],...
    'HorizontalAlignment','Right','Parent',htumorpanel);
htumornit = uicontrol('Style','text','String','',...
    'Position',[125,130,100,20],...
    'HorizontalAlignment','Left','Parent',htumorpanel);
%-------------------------------------------------------------------------%
htumorrec = uicontrol('Style','text','String','Recovered values',...
    'Position',[10,100,200,20],...
    'Parent',htumorpanel);
%-------------------------------------------------------------------------%
hkmftxtrec = uicontrol('Style','text','String','Kmf=',...
    'Position',[0,80,35,20],...
    'Parent',htumorpanel);
hkmfrec = uicontrol('Style','text','String','',...
    'Position',[35,80,70,20],...
    'Parent',htumorpanel);
%-------------------------------------------------------------------------%
hkfmtxtrec = uicontrol('Style','text','String','Kfm=',...
    'Position',[0,60,35,20],...
    'Parent',htumorpanel);
hkfmrec = uicontrol('Style','text','String','',...
    'Position',[35,60,70,20],...
    'Parent',htumorpanel);
%-------------------------------------------------------------------------%
hkbftxtrec = uicontrol('Style','text','String','Kbf=',...
    'Position',[110,80,35,20],...
    'Parent',htumorpanel);
hkbfrec = uicontrol('Style','text','String','',...
    'Position',[145,80,70,20],...
    'Parent',htumorpanel);
%-------------------------------------------------------------------------%
hkfbtxtrec = uicontrol('Style','text','String','Kfb=',...
    'Position',[110,60,35,20],...
    'Parent',htumorpanel);
hkfbrec = uicontrol('Style','text','String','',...
    'Position',[145,60,70,20],...
    'Parent',htumorpanel);
%-------------------------------------------------------------------------%
hrelerrtxt = uicontrol('Style','text','String','Relative error',...
    'Position',[5,30,210,20],...
    'Parent',htumorpanel);
hrelerr = uicontrol('Style','text','String','',...
    'Position',[5,10,210,20],...
    'Parent',htumorpanel);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
    function tumorigrand_Callback(hObject, evendata, handles)
        Kfbx=num2cell(1.5*rand(1,4));
        [Kfbx,Kmfx,Kbfx,Kfmx]=deal(Kfbx{:});
        set(hkfb,'String',num2str(Kfbx));
        set(hkmf,'String',num2str(Kmfx));
        set(hkbf,'String',num2str(Kbfx));
        set(hkfm,'String',num2str(Kfmx));
    end
%-------------------------------------------------------------------------%
    function tumorstartstop_Callback(hObject, evendata, handles)
        
        b=-1;
        startstop=~startstop;
        
        if startstop
            
            Kfbx=str2double(get(hkfb,'String'));
            Kmfx=str2double(get(hkmf,'String'));
            Kbfx=str2double(get(hkbf,'String'));
            Kfmx=str2double(get(hkfm,'String'));
            
            r=str2double(get(hreg,'String'));
            t = analysis_gui.DATA.ROI.Time_ROI./60;
            
            Ca = analysis_gui.IF.Averaged_IF;
            Cdata = analysis_gui.DATA.ROI.Averaged_ROI;            
            tbool=(t>t0);t=t(tbool);Ca=Ca(tbool);Cdata=Cdata(tbool);                      
            Ca=@(tt)(interp1([0;t],[0;Ca],tt,'linear',0)).';
            
            Ax=[[-(Kmfx+Kbfx);Kmfx],[Kfmx;-Kfmx]];
            Cx=concentration_K1(Ax,Ca,0,[0;0],t,glnodes,glweights);
            Cxdata=Kfbx*(alpha*Cx).';
                       
            nit=0;
            crit=Inf(5,5);
            relerrp = 0.5;
            
        end
        
        K_I= [Kfbx,Kmfx,Kbfx,Kfmx];
        
        while startstop&&any(crit(:)>critval)&&(nit<nitmax)%&&(relerrt>0.05)
            
            [Kfbx,Kmfx,Kbfx,Kfmx,Cx,Cxdata,relerr,nit,crit,relerrp]=...
                iterate_tumor_data(glnodes,glweights,t,alpha,Ca,Cdata,r,...
                Kfbx,Kmfx,Kbfx,Kfmx,Cx,Cxdata,relerr,nit,crit);
            
            set(hkfbrec,'String',num2str(Kfbx));
            set(hkmfrec,'String',num2str(Kmfx));
            set(hkbfrec,'String',num2str(Kbfx));
            set(hkfmrec,'String',num2str(Kfmx));                       
            set(hrelerr,'String',num2str(relerr));
            set(htumornit,'String',num2str(nit));
            
            
            while ( ~any( get(gcf,'Children')==ht ) )
                pause(1);
                fprintf('waiting...\n'); % Not necessary, but this helps to see if the code enters the loop
            end
            set(gcf,'CurrentAxes',ht);
            plot(t,Cdata,'r');
            hold on;
            plot(t,Cxdata,'g');
            xlabel('time [min]'); ylabel('concentration [kBq/cc]');
            title('Tumor');
            legend('Noisy data','Reconstructed data','Location','NorthEast');
            hold off;
            drawnow;
            
            b=0;
            
            if (nit>30)&&(relerrp<relerr)&&(relerrp<0.1)
                startstop =~startstop;
            end
            
            if (nit>30)&&(relerrp<relerr)&&(relerrp>0.1)
                startstop =~startstop;
                msgbox('Wrong initial guess. Choose other parameters');
                b=1;
            end
            
        end
        
        if (b==0)&&(relerrp<0.1)
            
            msgbox('Tumor converged');
            
            % rename and save variables
            K_parameters = {'Kfb','Kmf','Kbf','Kfm'};
            K = [Kfbx,Kmfx,Kbfx,Kfmx];
            Cxdata=Kfbx*(alpha*Cx).';
        
            DATA = analysis_gui.DATA;
            IF = analysis_gui.IF;
            KINETICS.parameters = K_parameters;
            KINETICS.values = K;
            FIT = struct('Fitting_curve',Cxdata);
            FIT.Relative_error = relerr;
            
            format short;
            aux_clock = fix(clock);
            date_time = strcat(num2str(aux_clock(1)),'-',num2str(aux_clock(2)),'-',num2str(aux_clock(3)),'_',num2str(aux_clock(4)),'_',num2str(aux_clock(5)),'_',num2str(aux_clock(6)));
            save([analysis_gui.OUTPUTfolder analysis_gui.slash char(analysis_gui.DATA.ROI.Study_Name) '_' analysis_gui.DATA.model '_' date_time '.mat'],'DATA','IF','KINETICS','FIT');
                       
        end
        
        if (b==0)&&(relerrp>0.1)
           msgbox('Wrong initial guess. Choose other parameters');
        end
        
        % If stop-button has not been pressed, the algorithm stops
        % because of the criterion --> change the startstopl here
        if startstop
            startstop=~startstop;
        end
               
    end

%-------------------------------------------------------------------------%
ht=axes('Units','Pixels','Position',[73.8 50 374 380]); 
%-------------------------------------------------------------------------%
set([h,...
    ht,htumorpanel,htumorig,...
    hkmftxt,hkmf,hkfmtxt,hkfm,hkbftxt,hkbf,hkfbtxt,hkfb,htumorigrand,...
    htumorstartstop,htumorrec,hkmftxtrec,hkmfrec,hkfmtxtrec,hkfmrec,...
    hkbftxtrec,hkbfrec,hkfbtxtrec,hkfbrec,htumornittxt,htumornit,...
    hrelerrtxt,hrelerr,hregtxttit,hregtxt,hreg],...
    'Units','normalized');
%-------------------------------------------------------------------------%
set(h,'Name','Tumor GUI');
set(h,'NumberTitle','off');
movegui(h,'center');
%-------------------------------------------------------------------------%
set(h,'Visible','on');
%-------------------------------------------------------------------------%

end