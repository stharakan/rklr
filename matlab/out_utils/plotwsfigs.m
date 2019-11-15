clear all
close all
datadir = './data/conv/';
files = dir(datadir);


for i = 1:length(files)
    % init, skip ahead if needed
    f = files(i);
    if f.isdir
        continue;
    end
    name = f.name;
    
    % inexact timings
    s = open([datadir,name]);
    ittimes = s.conv.ittimes;    
    ien_terrs = s.conv.terrs;
    ien_times = cumsum(ittimes);
    ien_leg = 'Inexact Newton';
    
    % exact newton
    enfile = [datadir,'en/',name];
    enflag = exist(enfile,'file');
    if enflag
        s = open(enfile);
        en_terrs = s.conv.terrs;
        ittimes = s.conv.ittimes;
        en_times = cumsum(ittimes);
        en_leg = 'Exact Newton';
    else
        en_terrs = [];
        en_times = [];
        en_leg = [];
    end
    
    % continuation
    wsfile = [datadir,'ws/',name];
    wsflag = exist(wsfile,'file');
    if wsflag
        s = open(wsfile);
        wsinfo = s.conv.wsinfo;
        iterrs = wsinfo.terrs;
        itimes = wsinfo.ittimes;
        if size(iterrs,2) ==2
            iterrs = iterrs(:,1);
            itimes = itimes(2:2:8);
        end
        
        fterrs = s.conv.terrs;
        ftimes = s.conv.ittimes;
        
        ws_only_times = cumsum(itimes);
        ws_only_times = [ws_only_times;ws_only_times(end)];
        ws_only_terrs = [iterrs;fterrs(1)];
        
        ittimes = [0;itimes;ftimes];
        ws_terrs = [ien_terrs(1);iterrs;fterrs];
        ws_times = cumsum(ittimes);
        ws_leg = 'm-Cont Inexact Newton';
    else
        ws_terrs = [];
        ws_times = [];
        ws_leg = [];
        ws_only_times = [];
        ws_only_terrs = [];
    end
    
    % steepest descent
    sdfile = [datadir,'sd/',name];
    sdflag = exist(sdfile,'file');
    maxtime = max([ws_times;en_times;ien_times]);
    if sdflag
        s = open(sdfile);
        sd_terrs = s.conv.terrs;
        ittimes = s.conv.ittimes;
        sd_times = cumsum(ittimes);
        sd_leg = 'Steepest Descent';
        
        % truncate if huge
        idx = sd_times<1.15 * maxtime;
        sd_terrs = sd_terrs(idx);
        sd_times = sd_times(idx);
    else
        sd_terrs = [];
        sd_times = [];
        sd_leg = [];
        
    end
    

    % continuation
    n2file = [datadir,'n2/',name];
    n2flag = exist(n2file,'file');
    n2flag = false;
    if n2flag
        s = open(n2file);
        wsinfo = s.convsm.wsinfo;
        iterrs = wsinfo.terrs;
        itimes = wsinfo.ittimes;
        if size(iterrs,2) ==2
            iterrs = iterrs(:,1);
            itimes = itimes(2:2:8);
        end
        
        fterrs = s.convsm.terrs;
        ftimes = s.convsm.ittimes;
        
        n2_only_times = cumsum(itimes);
        n2_only_times = [n2_only_times;n2_only_times(end)];
        n2_only_terrs = [iterrs;fterrs(1)];
        
        ittimes = [0;itimes;ftimes];
        n2_terrs = [ien_terrs(1);iterrs;fterrs];
        n2_times = cumsum(ittimes);
        n2_leg = 'm-Cont Inexact Newton at N/2';
    else
        n2_terrs = [];
        n2_times = [];
        n2_leg = [];
        n2_only_times = [];
        n2_only_terrs = [];
    end
    
    
    
    
    
    miny = min([min(ws_terrs);min(en_terrs);min(ien_terrs)]);
    maxy = max([max(ws_terrs);max(en_terrs);max(ien_terrs)]);
    figure;
    semilogy(sd_times,sd_terrs,en_times,en_terrs,ien_times,...
        ien_terrs,ws_times,ws_terrs,n2_times,n2_terrs,'LineWidth',2);
    hold on
    semilogy(ws_only_times,ws_only_terrs,'ko','Linewidth',5);
    grid on
    if sdflag
        if enflag
            if wsflag
                if n2flag
                    legend(sd_leg,en_leg,ien_leg,ws_leg,n2_leg);
                else
                    legend(sd_leg,en_leg,ien_leg,ws_leg);
                end
            else
                if n2flag
                    legend(sd_leg,en_leg,ien_leg,n2_leg);
                else
                    legend(sd_leg,en_leg,ien_leg);
                end
            end
        else
            if wsflag
                if n2flag
                    legend(sd_leg,ien_leg,ws_leg,n2_leg);
                else
                    legend(sd_leg,ien_leg,ws_leg);
                end
            else
                if n2flag
                    legend(sd_leg,ien_leg,n2_leg);
                else
                    legend(sd_leg,ien_leg);
                end
            end
        end
    else
        if enflag
            if wsflag
                if n2flag
                    legend(en_leg,ien_leg,ws_leg,n2_leg);
                else
                    legend(en_leg,ien_leg,ws_leg);
                end
            else
                if n2flag
                    legend(en_leg,ien_leg,n2_leg);
                else
                    legend(en_leg,ien_leg);
                end
            end
        else
            if wsflag
                if n2flag
                    legend(ien_leg,ws_leg,n2_leg);
                else
                    legend(ien_leg,ws_leg);
                end
            else
                if n2flag
                    legend(ien_leg,n2_leg);
                else
                    legend(ien_leg);
                end
            end
        end
    end
    
    
    axis( [0 maxtime*1.15 miny*0.8 maxy * 1.15]);
    axis 'auto x'
    
    xlabel('Time(s)','FontSize',20);
    ylabel('Testing error','FontSize',20);
    title([name,' Error vs time'],'FontSize',20);
    hold off
        
end