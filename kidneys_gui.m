function kidneys_gui

global analysis_gui
%-------------------------------------------------------------------------%
% whitecolor=[1,1,1];
startstop=false;
critval=1e-4;
relerr=Inf;
nitmax=500;
%Vb=; Vt=;
[glnodes,glweights]=gauss_legendre(4);
alpha=[1,1,1];
%-------------------------------------------------------------------------%
t0=0;
%-------------------------------------------------------------------------%
[t,Ca,Cdata]=deal(0);
[Kfax,Kmax,Kafx,Kmfx,Kfmx,Ktmx,Kutx,Cx1,Cx2,Cx]=deal(0);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
h = figure('Visible','off','Position',[200,200,800,450],...
    'MenuBar','none','Toolbar','figure');
%-------------------------------------------------------------------------%

% KIDNEYS
%-------------------------------------------------------------------------%
hkidpanel = uipanel('Title','Kidneys',...
    'TitlePosition','centertop',...
    'Position',[460/800,63/450,330/800,360/450]);
%-------------------------------------------------------------------------%
hkidig = uicontrol('Style','text','String','Initial guess',...
    'Position',[5,315,110,20],...
    'HorizontalAlignment','Left','Parent',hkidpanel);
%-------------------------------------------------------------------------%
hkidigrand = uicontrol('Style','pushbutton','String','Random',...
    'Position',[220,315,105,25],...
    'HorizontalAlignment','Right','Parent',hkidpanel,...
    'Callback',{@kidigrand_Callback});
%-------------------------------------------------------------------------%
hktmtxt = uicontrol('Style','text','String','Ktm=',...
    'Position',[0,280,35,20],...
    'Parent',hkidpanel);
hktm = uicontrol('Style','edit','String','',...
    'Position',[35,280,70,25],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hkuttxt = uicontrol('Style','text','String','Kut=',...
    'Position',[0,250,35,20],...
    'Parent',hkidpanel);
hkut = uicontrol('Style','edit','String','',...
    'Position',[35,250,70,25],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hkmftxt = uicontrol('Style','text','String','Kmf=',...
    'Position',[110,280,35,20],...
    'Parent',hkidpanel);
hkmf = uicontrol('Style','edit','String','',...
    'Position',[145,280,70,25],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hkfmtxt = uicontrol('Style','text','String','Kfm=',...
    'Position',[110,250,35,20],...
    'Parent',hkidpanel);
hkfm = uicontrol('Style','edit','String','',...
    'Position',[145,250,70,25],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hkaftxt = uicontrol('Style','text','String','Kaf=',...
    'Position',[220,280,35,20],...
    'Parent',hkidpanel);
hkaf = uicontrol('Style','edit','String','',...
    'Position',[255,280,70,25],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hkfatxt = uicontrol('Style','text','String','Kfa=',...
    'Position',[220,250,35,20],...
    'Parent',hkidpanel);
hkfa = uicontrol('Style','edit','String','',...
    'Position',[255,250,70,25],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hkmatxt = uicontrol('Style','text','String','Kma=',...
    'Position',[220,220,35,20],...
    'Parent',hkidpanel);
hkma = uicontrol('Style','edit','String','',...
    'Position',[255,220,70,25],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hregtxttit = uicontrol('Style','text','String','Regularization parameter',...
    'Position',[5,215,210,20],...
    'Parent',hkidpanel);
hregtxt = uicontrol('Style','text','String','r=',...
    'Position',[70,190,30,20],...
    'Parent',hkidpanel);
hreg = uicontrol('Style','edit','String','1e4',...
    'Position',[100,190,70,25],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hkidstartstop = uicontrol('Style','pushbutton','String',...
    'Start/Stop iterations','Position',[10,155,200,25],...
    'Parent',hkidpanel,...
    'Callback',{@kidstartstop_Callback});
%-------------------------------------------------------------------------%
hkidnittxt = uicontrol('Style','text','String','Iteration',...
    'Position',[15,130,105,20],...
    'HorizontalAlignment','Right','Parent',hkidpanel);
hkidnit = uicontrol('Style','text','String','',...
    'Position',[125,130,100,20],...
    'HorizontalAlignment','Left','Parent',hkidpanel);
%-------------------------------------------------------------------------%
hkidrec = uicontrol('Style','text','String','Recovered values',...
    'Position',[10,100,200,20],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hktmtxtrec = uicontrol('Style','text','String','Ktm=',...
    'Position',[0,80,35,20],...
    'Parent',hkidpanel);
hktmrec = uicontrol('Style','text','String','',...
    'Position',[35,80,70,20],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hkuttxtrec = uicontrol('Style','text','String','Kut=',...
    'Position',[0,60,35,20],...
    'Parent',hkidpanel);
hkutrec = uicontrol('Style','text','String','',...
    'Position',[35,60,70,20],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hkmftxtrec = uicontrol('Style','text','String','Kmf=',...
    'Position',[110,80,35,20],...
    'Parent',hkidpanel);
hkmfrec = uicontrol('Style','text','String','',...
    'Position',[145,80,70,20],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hkfmtxtrec = uicontrol('Style','text','String','Kfm=',...
    'Position',[110,60,35,20],...
    'Parent',hkidpanel);
hkfmrec = uicontrol('Style','text','String','',...
    'Position',[145,60,70,20],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hkaftxtrec = uicontrol('Style','text','String','Kaf=',...
    'Position',[220,80,35,20],...
    'Parent',hkidpanel);
hkafrec = uicontrol('Style','text','String','',...
    'Position',[255,80,70,20],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hkfatxtrec = uicontrol('Style','text','String','Kfa=',...
    'Position',[220,60,35,20],...
    'Parent',hkidpanel);
hkfarec = uicontrol('Style','text','String','',...
    'Position',[255,60,70,20],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hkmatxtrec = uicontrol('Style','text','String','Kma=',...
    'Position',[220,40,35,20],...
    'Parent',hkidpanel);
hkmarec = uicontrol('Style','text','String','',...
    'Position',[255,40,70,20],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%
hrelerrtxt = uicontrol('Style','text','String','Relative error',...
    'Position',[5,30,210,20],...
    'Parent',hkidpanel);
hrelerr = uicontrol('Style','text','String','',...
    'Position',[5,10,210,20],...
    'Parent',hkidpanel);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
    function kidigrand_Callback(hObject, evendata, handles)
        Kfax=num2cell(1.5*rand(1,7));
        [Kfax,Kmax,Kafx,Kmfx,Kfmx,Ktmx,Kutx]=deal(Kfax{:});
        set(hkfa,'String',num2str(Kfax));
        set(hkma,'String',num2str(Kmax));
        set(hkaf,'String',num2str(Kafx));
        set(hkmf,'String',num2str(Kmfx));
        set(hkfm,'String',num2str(Kfmx));
        set(hktm,'String',num2str(Ktmx));     
        set(hkut,'String',num2str(Kutx));
    end
%-------------------------------------------------------------------------%
    function kidstartstop_Callback(hObject, evendata, handles)
        
        b=-1;
        startstop=~startstop;
        
        if startstop
            
            Kfax=str2double(get(hkfa,'String'));
            Kmax=str2double(get(hkma,'String'));
            Kafx=str2double(get(hkaf,'String'));
            Kmfx=str2double(get(hkmf,'String'));
            Kfmx=str2double(get(hkfm,'String'));
            Ktmx=str2double(get(hktm,'String'));
            Kutx=str2double(get(hkut,'String'));
            
            r=str2double(get(hreg,'String'));
            t = analysis_gui.DATA.ROI.Time_ROI./60;
            Ca = analysis_gui.IF.Averaged_IF;
            Cdata = analysis_gui.DATA.ROI.Averaged_ROI;
            tbool=(t>t0);t=t(tbool);Ca=Ca(tbool);Cdata=Cdata(tbool);
            Ca=@(tt)(interp1([0;t],[0;Ca],tt,'linear',0)).';
            
            Ax=[[-(Kafx+Kmfx);Kmfx;0],[Kfmx;-(Kfmx+Ktmx);Ktmx],[0;0;-Kutx]];
            Cx1=concentration_K1_kid(Ax,Ca,[1;0;0],0,[0;0;0],t,glnodes,glweights); %[Ca;0;0] on e1
            Cx2=concentration_K1_kid(Ax,Ca,[0;1;0],0,[0;0;0],t,glnodes,glweights); %[0;Ca;0] on e2
            Cx=Kfax*Cx1+Kmax*Cx2;
            Cxdata=(alpha*Cx).';
            
            nit=0;
            crit=Inf(8,8);
            relerrp = 0.5;
            
        end
        
        K_I = [Kfax,Kmax,Kafx,Kmfx,Kfmx,Ktmx,Kutx];
        
        while startstop&&any(crit(:)>critval)&&(nit<nitmax)%&&(relerrt>0.05)
            
            [Kfax,Kmax,Kafx,Kmfx,Kfmx,Ktmx,Kutx,Cx1,Cx2,Cxdata,relerr,nit,crit,relerrp]=...
                iterate_kid_data(glnodes,glweights,t,alpha,Ca,Cdata,r,...
                Kfax,Kmax,Kafx,Kmfx,Kfmx,Ktmx,Kutx,Cx1,Cx2,Cxdata,relerr,nit,crit);
                        
            set(hkfarec,'String',num2str(Kfax)); 
            set(hkmarec,'String',num2str(Kmax));
            set(hkafrec,'String',num2str(Kafx));
            set(hkmfrec,'String',num2str(Kmfx));
            set(hkfmrec,'String',num2str(Kfmx));
            set(hktmrec,'String',num2str(Ktmx));
            set(hkutrec,'String',num2str(Kutx));
            set(hrelerr,'String',num2str(relerr));
            set(hkidnit,'String',num2str(nit));
            
            set(gcf,'CurrentAxes',hk);
            plot(t,Cdata,'r');
            hold on;
            plot(t,Cxdata,'g');
            xlabel('time [min]'); ylabel('concentration [kBq/cc]');
            title('Kidneys');
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
            
            msgbox('Kidneys converged');
            
            % rename and save variables
            K_parameters = {'Kfa','Kma','Kaf','Kmf','Kfm','Ktm','Kut'};
            K = [Kfax,Kmax,Kafx,Kmfx,Kfmx,Ktmx,Kutx];
            Cx=Kfax*Cx1+Kmax*Cx2;
            Cxdata=(alpha*Cx).';
            
            DATA = analysis_gui.DATA;
            IF = analysis_gui.IF;
            KINETICS.parameters = K_parameters;
            KINETICS.values = K;
            FIT.Fitting_curve = Cxdata;
            FIT.Relative_error = relerr;
            
            format short;
            aux_clock = fix(clock);
            date_time = strcat(num2str(aux_clock(1)),'-',num2str(aux_clock(2)),'-',num2str(aux_clock(3)),'_',num2str(aux_clock(4)),'_',num2str(aux_clock(5)),'_',num2str(aux_clock(6)));
            save([analysis_gui.OUTPUTfolder analysis_gui.slash char(analysis_gui.DATA.ROI.Study_Name) '_' analysis_gui.DATA.model '_' date_time '.mat'],'DATA','IF','KINETICS','FIT');
                        
        end
        
        if (b==0)&&(relerrp>0.1)
           msgbox('Wrong initial guess. Choose other parameters');
        end

        % If stop has not been pressed, that is, the algorithm stopped
        % because of the criterion, we change the startstopl here
        if startstop
            startstop=~startstop;
        end
        
        
    end

%-------------------------------------------------------------------------%
hk=axes('Units','Pixels','Position',[73.8 50 374 380]); 
%-------------------------------------------------------------------------%
set([h,...
    hk,hkidpanel,hkidig,...
    hktmtxt,hktm,hkuttxt,hkut,hkmftxt,hkmf,hkfmtxt,hkfm,hkaftxt,hkaf,hkfatxt,hkfa,hkmatxt,hkma,...
    hkidigrand,hkidstartstop,hkidrec,...
    hktmtxtrec,hktmrec,hkuttxtrec,hkutrec,hkmftxtrec,hkmfrec,hkfmtxtrec,hkfmrec,...
    hkaftxtrec,hkafrec,hkfatxtrec,hkfarec,hkmatxtrec,hkmarec,hkmatxtrec,hkmarec,...
    hkidnittxt,hkidnit,hrelerrtxt,hrelerr,hregtxttit,hregtxt,hreg],...
    'Units','normalized');
%-------------------------------------------------------------------------%
set(h,'Name','Kidneys GUI');
set(h,'NumberTitle','off');
movegui(h,'center');
%-------------------------------------------------------------------------%
set(h,'Visible','on');
%-------------------------------------------------------------------------%

end