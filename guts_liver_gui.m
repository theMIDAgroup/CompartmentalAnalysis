function guts_liver_gui

global analysis_gui
%-------------------------------------------------------------------------%
% whitecolor=[1,1,1];
startstopg=false;
critvalg=1e-4;
relerrg=Inf;
nitgmax=500;
startstopl=false;
critvall=1e-4;
relerrl=Inf;
nitlmax=500;
V=0.3; %Vb blood fraction
[glnodes,glweights]=gauss_legendre(4);
alpha=[1,1];
%-------------------------------------------------------------------------%
t0=0;
%-------------------------------------------------------------------------%
[t,Ca,Cdatag,Cdatal]=deal(0);
[Kgax,Ktgx,Kfpx,Kgtx,Cxg]=deal(0);
[Kfax,Kmfx,Ksfx,Kfmx,Cxl]=deal(0);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
h = figure('Visible','off','Position',[200,200,980,500],...
    'MenuBar','none','Toolbar','figure');
%-------------------------------------------------------------------------%

% GUT
%-------------------------------------------------------------------------%
hgutspanel = uipanel('Title','Gut',...
    'TitlePosition','centertop',...
    'Position',[520/980,83/500,220/980,360/500]);
%-------------------------------------------------------------------------%
hgutsig = uicontrol('Style','text','String','Initial guess',...
    'Position',[5,315,110,20],...
    'HorizontalAlignment','Left','Parent',hgutspanel);
%-------------------------------------------------------------------------%
hgutsigrand = uicontrol('Style','pushbutton','String','Random',...
    'Position',[110,315,105,25],...
    'HorizontalAlignment','Right','Parent',hgutspanel,...
    'Callback',{@gutsigrand_Callback});
%-------------------------------------------------------------------------%
hktgtxt = uicontrol('Style','text','String','Ktg=',...
    'Position',[0,280,35,20],...
    'Parent',hgutspanel);
hktg = uicontrol('Style','edit','String','',...
    'Position',[35,280,70,25],...
    'Parent',hgutspanel);
%-------------------------------------------------------------------------%
hkgttxt = uicontrol('Style','text','String','Kgt=',...
    'Position',[0,250,35,20],...
    'Parent',hgutspanel);
hkgt = uicontrol('Style','edit','String','',...
    'Position',[35,250,70,25],...
    'Parent',hgutspanel);
%-------------------------------------------------------------------------%
hkvgtxt = uicontrol('Style','text','String','Kfp=',...
    'Position',[110,280,35,20],...
    'Parent',hgutspanel);
hkvg = uicontrol('Style','edit','String','',...
    'Position',[145,280,70,25],...
    'Parent',hgutspanel);
%-------------------------------------------------------------------------%
hkgatxt = uicontrol('Style','text','String','Kga=',...
    'Position',[110,250,35,20],...
    'Parent',hgutspanel);
hkga = uicontrol('Style','edit','String','',...
    'Position',[145,250,70,25],...
    'Parent',hgutspanel);
%-------------------------------------------------------------------------%
hreggtxttit = uicontrol('Style','text','String','Regularization parameter',...
    'Position',[5,215,210,20],...
    'Parent',hgutspanel);
hreggtxt = uicontrol('Style','text','String','r=',...
    'Position',[70,190,30,20],...
    'Parent',hgutspanel);
hregg = uicontrol('Style','edit','String','1e4',...
    'Position',[100,190,70,25],...
    'Parent',hgutspanel);
%-------------------------------------------------------------------------%
hgutsstartstop = uicontrol('Style','pushbutton','String',...
    'Start/Stop iterations','Position',[10,155,200,25],...
    'Parent',hgutspanel,...
    'Callback',{@gutsstartstop_Callback});
%-------------------------------------------------------------------------%
hgutsnittxt = uicontrol('Style','text','String','Iteration',...
    'Position',[15,130,105,20],...
    'HorizontalAlignment','Right','Parent',hgutspanel);
hgutsnit = uicontrol('Style','text','String','',...
    'Position',[125,130,100,20],...
    'HorizontalAlignment','Left','Parent',hgutspanel);
%-------------------------------------------------------------------------%
hgutsrec = uicontrol('Style','text','String','Recovered values',...
    'Position',[10,100,200,20],...
    'Parent',hgutspanel);
%-------------------------------------------------------------------------%
hktgtxtrec = uicontrol('Style','text','String','Ktg=',...
    'Position',[0,80,35,20],...
    'Parent',hgutspanel);
hktgrec = uicontrol('Style','text','String','',...
    'Position',[35,80,70,20],...
    'Parent',hgutspanel);
%-------------------------------------------------------------------------%
hkgttxtrec = uicontrol('Style','text','String','Kgt=',...
    'Position',[0,60,35,20],...
    'Parent',hgutspanel);
hkgtrec = uicontrol('Style','text','String','',...
    'Position',[35,60,70,20],...
    'Parent',hgutspanel);
%-------------------------------------------------------------------------%
hkvgtxtrec = uicontrol('Style','text','String','Kfp=',...
    'Position',[110,80,35,20],...
    'Parent',hgutspanel);
hkvgrec = uicontrol('Style','text','String','',...
    'Position',[145,80,70,20],...
    'Parent',hgutspanel);
%-------------------------------------------------------------------------%
hkgatxtrec = uicontrol('Style','text','String','Kga=',...
    'Position',[110,60,35,20],...
    'Parent',hgutspanel);
hkgarec = uicontrol('Style','text','String','',...
    'Position',[145,60,70,20],...
    'Parent',hgutspanel);
%-------------------------------------------------------------------------%
hrelerrgtxt = uicontrol('Style','text','String','Relative error',...
    'Position',[5,30,210,20],...
    'Parent',hgutspanel);
hrelerrg = uicontrol('Style','text','String','',...
    'Position',[5,10,210,20],...
    'Parent',hgutspanel);
%-------------------------------------------------------------------------%

%LIVER
%-------------------------------------------------------------------------%
hliverpanel = uipanel('Title','Liver',...
    'TitlePosition','centertop',...
    'Position',[750/980,83/500,220/980,360/500]);
%-------------------------------------------------------------------------%
hliverig = uicontrol('Style','text','String','Initial guess',...
    'Position',[5,315,110,20],...
    'HorizontalAlignment','Left','Parent',hliverpanel);
%-------------------------------------------------------------------------%
hliverigrand = uicontrol('Style','pushbutton','String','Random',...
    'Position',[110,315,105,25],...
    'HorizontalAlignment','Right','Parent',hliverpanel,...
    'Callback',{@liverigrand_Callback});
%-------------------------------------------------------------------------%
hkmftxt = uicontrol('Style','text','String','Kmf=',...
    'Position',[0,280,35,20],...
    'Parent',hliverpanel);
hkmf = uicontrol('Style','edit','String','',...
    'Position',[35,280,70,25],...
    'Parent',hliverpanel);
%-------------------------------------------------------------------------%
hkfmtxt = uicontrol('Style','text','String','Kfm=',...
    'Position',[0,250,35,20],...
    'Parent',hliverpanel);
hkfm = uicontrol('Style','edit','String','',...
    'Position',[35,250,70,25],...
    'Parent',hliverpanel);
%-------------------------------------------------------------------------%
hkpftxt = uicontrol('Style','text','String','Ksf=',...
    'Position',[110,280,35,20],...
    'Parent',hliverpanel);
hkpf = uicontrol('Style','edit','String','',...
    'Position',[145,280,70,25],...
    'Parent',hliverpanel);
%-------------------------------------------------------------------------%
hklatxt = uicontrol('Style','text','String','Kfa=',...
    'Position',[110,250,35,20],...
    'Parent',hliverpanel);
hkla = uicontrol('Style','edit','String','',...
    'Position',[145,250,70,25],...
    'Parent',hliverpanel);
%-------------------------------------------------------------------------%
hregltxttit = uicontrol('Style','text','String','Regularization parameter',...
    'Position',[5,215,210,20],...
    'Parent',hliverpanel);
hregltxt = uicontrol('Style','text','String','r=',...
    'Position',[70,190,30,20],...
    'Parent',hliverpanel);
hregl = uicontrol('Style','edit','String','1e4',...
    'Position',[100,190,70,25],...
    'Parent',hliverpanel);
%-------------------------------------------------------------------------%
hliverstartstop = uicontrol('Style','pushbutton','String',...
    'Start/Stop iterations','Position',[10,155,200,25],...
    'Parent',hliverpanel,...
    'Callback',{@liverstartstop_Callback});
%-------------------------------------------------------------------------%
hlivernittxt = uicontrol('Style','text','String','Iteration',...
    'Position',[15,130,105,20],...
    'HorizontalAlignment','Right','Parent',hliverpanel);
hlivernit = uicontrol('Style','text','String','',...
    'Position',[125,130,100,20],...
    'HorizontalAlignment','Left','Parent',hliverpanel);
%-------------------------------------------------------------------------%
hliverrec = uicontrol('Style','text','String','Recovered values',...
    'Position',[10,100,200,20],...
    'Parent',hliverpanel);
%-------------------------------------------------------------------------%
hkmftxtrec = uicontrol('Style','text','String','Kmf=',...
    'Position',[0,80,35,20],...
    'Parent',hliverpanel);
hkmfrec = uicontrol('Style','text','String','',...
    'Position',[35,80,70,20],...
    'Parent',hliverpanel);
%-------------------------------------------------------------------------%
hkfmtxtrec = uicontrol('Style','text','String','Kfm=',...
    'Position',[0,60,35,20],...
    'Parent',hliverpanel);
hkfmrec = uicontrol('Style','text','String','',...
    'Position',[35,60,70,20],...
    'Parent',hliverpanel);
%-------------------------------------------------------------------------%
hkpftxtrec = uicontrol('Style','text','String','Ksf=',...
    'Position',[110,80,35,20],...
    'Parent',hliverpanel);
hkpfrec = uicontrol('Style','text','String','',...
    'Position',[145,80,70,20],...
    'Parent',hliverpanel);
%-------------------------------------------------------------------------%
hklatxtrec = uicontrol('Style','text','String','Kfa=',...
    'Position',[110,60,35,20],...
    'Parent',hliverpanel);
hklarec = uicontrol('Style','text','String','',...
    'Position',[145,60,70,20],...
    'Parent',hliverpanel);
%-------------------------------------------------------------------------%
hrelerrltxt = uicontrol('Style','text','String','Relative error',...
    'Position',[5,30,210,20],...
    'Parent',hliverpanel);
hrelerrl = uicontrol('Style','text','String','',...
    'Position',[5,10,210,20],...
    'Parent',hliverpanel);

%-------------------------------------------------------------------------%
    function gutsigrand_Callback(hObject, evendata, handles)
        Kgax=num2cell(1.5*rand(1,4));
        [Kgax,Ktgx,Kfpx,Kgtx]=deal(Kgax{:});
        set(hkga,'String',num2str(Kgax));
        set(hktg,'String',num2str(Ktgx));
        set(hkvg,'String',num2str(Kfpx));
        set(hkgt,'String',num2str(Kgtx));
    end
%-------------------------------------------------------------------------%
    function gutsstartstop_Callback(hObject, evendata, handles)
        
        b_g=-1;
        startstopg=~startstopg;
        
        if startstopg
            
            Kgax=str2double(get(hkga,'String'));
            Ktgx=str2double(get(hktg,'String'));
            Kfpx=str2double(get(hkvg,'String'));
            Kgtx=str2double(get(hkgt,'String'));
            
            rg=str2double(get(hregg,'String'));
            t = analysis_gui.DATA.ROI.Time_ROI./60;
            
            Ca = analysis_gui.IF.Averaged_IF;
            Cdatag = analysis_gui.GUT.Averaged_GUT;
            Cdatal = analysis_gui.DATA.ROI.Averaged_ROI;
            tbool=(t>t0);t=t(tbool);Ca=Ca(tbool);Cdatag=Cdatag(tbool);Cdatal=Cdatal(tbool);
            Ca=@(tt)(interp1([0;t],[0;Ca],tt,'linear',0)).';
            
            Agx=[[-(Ktgx+Kfpx);Ktgx],[Kgtx;-Kgtx]];
            Cxg=concentration_K1(Agx,Ca,0,[0;0],t,glnodes,glweights);
            Cxdatag=Kgax*(alpha*Cxg).';
            
            nitg=0;
            critg=Inf(5,5);
            relerrpg = 0.5;
            
        end
        
        Kg_I= [Kgax,Ktgx,Kfpx,Kgtx];
        
        while startstopg&&any(critg(:)>critvalg)&&(nitg<nitgmax)%&&(relerrpg>0.05)
            
            [Kgax,Ktgx,Kfpx,Kgtx,Cxg,Cxdatag,relerrg,nitg,critg, relerrpg]=...
                iterate_guts_data(glnodes,glweights,t,alpha,Ca,Cdatag,rg,...
                Kgax,Ktgx,Kfpx,Kgtx,Cxg,Cxdatag,relerrg,nitg,critg);
            
            set(hkgarec,'String',num2str(Kgax));
            set(hktgrec,'String',num2str(Ktgx));
            set(hkvgrec,'String',num2str(Kfpx));
            set(hkgtrec,'String',num2str(Kgtx));
            set(hrelerrg,'String',num2str(relerrg));
            set(hgutsnit,'String',num2str(nitg));
            
            set(gcf,'CurrentAxes',hg);
            plot(t,Cdatag,'r');
            hold on;
            plot(t,Cxdatag,'g');
            xlabel('time [min]'); ylabel('concentration [kBq/cc]');
            title('Gut');
            legend('Noisy data','Reconstructed data','Location','NorthEast');
            hold off;
            drawnow;
            
            b_g=0;
            
            if (nitg>30)&&(relerrpg<relerrg)&&(relerrpg<0.1)
                startstopg =~startstopg;
            end
            
            if (nitg>30)&&(relerrpg<relerrg)&&(relerrpg>0.1)
                startstopg =~startstopg;
                msgbox('Wrong initial guess. Choose other parameters');
                b_g=1;
            end
            
        end
        
        if (b_g==0)&&(relerrpg<0.1)
            
            msgbox('Gut converged');
            
            % rename and save variables
            Kg_parameters = {'Kga','Ktg','Kfp','Kgt'};
            Kg = [Kgax,Ktgx,Kfpx,Kgtx];
            Cxdatag=Kgax*(alpha*Cxg).';
            
            DATA = analysis_gui.DATA;
            GUT = analysis_gui.GUT;
            IF = analysis_gui.IF;
            KINETICS_GUT.parameters = Kg_parameters;
            KINETICS_GUT.values = Kg;
            FIT_GUT = struct('Fitting_curve',Cxdatag);
            FIT_GUT.Relative_error = relerrg;
            
            format short;
            aux_clock = fix(clock);
            date_time = strcat(num2str(aux_clock(1)),'-',num2str(aux_clock(2)),'-',num2str(aux_clock(3)),'_',num2str(aux_clock(4)),'_',num2str(aux_clock(5)),'_',num2str(aux_clock(6)));
            save([analysis_gui.OUTPUTfolder analysis_gui.slash char(analysis_gui.DATA.ROI.Study_Name) '_Gut_' date_time '.mat'],'DATA','GUT','IF','KINETICS_GUT','FIT_GUT');
            
        end
        
        if (b_g==0)&&(relerrpg>0.1)
           msgbox('Wrong initial guess. Choose other parameters');
        end
               
        % If stop has not been pressed, that is, the algorithm stopped
        % because of the criterion, we change the startstopl here
        if startstopg
            startstopg=~startstopg;
        end
        
               
    end

%-------------------------------------------------------------------------%
    function  liverigrand_Callback(hObject, evendata, handles)
        Kfax=num2cell(1.5*rand(1,4));
        [Kfax,Kmfx,Ksfx,Kfmx]=deal(Kfax{:});
        set(hkla,'String',num2str(Kfax));
        set(hkmf,'String',num2str(Kmfx));
        set(hkpf,'String',num2str(Ksfx));
        set(hkfm,'String',num2str(Kfmx));        
    end
%-------------------------------------------------------------------------%
    function liverstartstop_Callback(hObject, evendata, handles)
        
        b=-1;
        startstopl=~startstopl;
        
        if startstopl
            
            Kfax=str2double(get(hkla,'String'));
            Kmfx=str2double(get(hkmf,'String'));
            Ksfx=str2double(get(hkpf,'String'));
            Kfmx=str2double(get(hkfm,'String'));
            
            Alx=[[-(Kmfx+Ksfx);Kmfx],[Kfmx;-Kfmx]];
            rl=str2double(get(hregl,'String'));
            Agx=[[-(Ktgx+Kfpx);Ktgx],[Kgtx;-Kgtx]];
            Cv=concentration_Cv(Agx,Ca,t,Cxg,glnodes,glweights).';
            Cv=@(tt)(interp1([0;t],[0;Cv],tt,'linear',0)).';
            Cxla=concentration_K1(Alx,Ca,0,[0;0],t,glnodes,glweights);
            Cxlv=concentration_K1(Alx,Cv,0,[0;0],t,glnodes,glweights);
            Cxl=Kfax*Cxla+Kfpx*Cxlv;
            Cxdatal=(1-V)*(alpha*Cxl).'+V/4*(Ca(t)+3*Cv(t)).';
            
            nitl=0;
            critl=Inf(5,5);
            
        end
        
        Kl_I = [Kfax,Kmfx,Ksfx,Kfmx];
        
        while startstopl&&any(critl(:)>critvall)&&(nitl<nitlmax)
            
            [Kfax,Kmfx,Ksfx,Kfmx,Cxla,Cxlv,Cxdatal,relerrl,nitl,critl,relerrpl]=...
                iterate_liver_data(glnodes,glweights,t,alpha,V,Ca,Kfpx,Cv,Cdatal,rl,...
                Kfax,Kmfx,Ksfx,Kfmx,Cxla,Cxlv,Cxdatal,relerrl,nitl,critl);
            
            set(hklarec,'String',num2str(Kfax));
            set(hkmfrec,'String',num2str(Kmfx));
            set(hkpfrec,'String',num2str(Ksfx));
            set(hkfmrec,'String',num2str(Kfmx));
            set(hrelerrl,'String',num2str(relerrl));
            set(hlivernit,'String',num2str(nitl));
            
            set(gcf,'CurrentAxes',hl);
            plot(t,Cdatal,'r');
            hold on;
            plot(t,Cxdatal,'g');
            xlabel('time [min]'); ylabel('concentration [kBq/cc]');
            title('Liver');
            legend('Noisy data','Reconstructed data','Location','NorthEast');
            hold off;
            drawnow;
            
            b=0;
            
            if (nitl>30)&&(relerrpl<relerrl)&&(relerrpl<0.1)
                startstopl =~startstopl;
            end
            
            if (nitl>30)&&(relerrpl<relerrl)&&(relerrpl>0.1)
                startstopl =~startstopl;
                msgbox('Wrong initial guess. Choose other parameters');
                b=1;
            end
            
        end
        
        if (b==0)&&(relerrpl<0.1)
            
            msgbox('Liver converged');
            
            % rename and save variables
            Kl_parameters = {'Kfa','Kmf','Ksf','Kfm'};
            Kl = [Kfax,Kmfx,Ksfx,Kfmx];
            
            Agx=[[-(Ktgx+Kfpx);Ktgx],[Kgtx;-Kgtx]];
            Cv=concentration_Cv(Agx,Ca,t,Cxg,glnodes,glweights).';
            Cv=@(tt)(interp1([0;t],[0;Cv],tt,'linear',0)).';
            Cxdatal=(1-V)*(alpha*Cxl).'+V/4*(Ca(t)+3*Cv(t)).';
            
            DATA = analysis_gui.DATA;
            GUT = analysis_gui.GUT;
            IF = analysis_gui.IF;
            KINETICS_LIVER.parameters = Kl_parameters;
            KINETICS_LIVER.values = Kl;
            FIT_LIVER = struct('Fitting_curve',Cxdatal);
            FIT_LIVER.Relative_error = relerrl;
            
            format short;
            aux_clock = fix(clock);
            date_time = strcat(num2str(aux_clock(1)),'-',num2str(aux_clock(2)),'-',num2str(aux_clock(3)),'_',num2str(aux_clock(4)),'_',num2str(aux_clock(5)),'_',num2str(aux_clock(6)));
            save([analysis_gui.OUTPUTfolder analysis_gui.slash char(analysis_gui.DATA.ROI.Study_Name) '_' analysis_gui.DATA.model '_' date_time '.mat'],'DATA','GUT','IF','KINETICS_LIVER','FIT_LIVER');
            
        end
        
         if (b==0)&&(relerrpl>0.1)
           msgbox('Wrong initial guess. Choose other parameters');
         end
        
        % If stop has not been pressed, that is, the algorithm stopped
        % because of the criterion, we change the startstopl here
        if startstopl
            startstopl=~startstopl;
        end
                 
    end

%-------------------------------------------------------------------------%
hg=axes('Units','Pixels','Position',[73.8 300 434 171.15]);
hl=axes('Units','Pixels','Position',[73.8 48 434 171.15]);
%-------------------------------------------------------------------------%
set([h,...
    hg,hgutspanel,hgutsig,...
    hktgtxt,hktg,hkgttxt,hkgt,hkvgtxt,hkvg,hkgatxt,hkga,hgutsigrand,...
    hgutsstartstop,hgutsrec,hktgtxtrec,hktgrec,hkgttxtrec,hkgtrec,...
    hkvgtxtrec,hkvgrec,hkgatxtrec,hkgarec,hgutsnittxt,hgutsnit,...
    hrelerrgtxt,hrelerrg,hreggtxttit,hreggtxt,hregg,...
    hl,hliverpanel,hliverig,...
    hkmftxt,hkmf,hkfmtxt,hkfm,hkpftxt,hkpf,hklatxt,hkla,hliverigrand,...
    hliverstartstop,hliverrec,hkmftxtrec,hkmfrec,hkfmtxtrec,hkfmrec,...
    hkpftxtrec,hkpfrec,hklatxtrec,hklarec,hlivernittxt,hlivernit,...
    hrelerrltxt,hrelerrl,hregltxttit,hregltxt,hregl],...
    'Units','normalized');
%-------------------------------------------------------------------------%
set(h,'Name','Gut-Liver GUI');
set(h,'NumberTitle','off');
movegui(h,'center');
%-------------------------------------------------------------------------%
set(h,'Visible','on');
%-------------------------------------------------------------------------%

end