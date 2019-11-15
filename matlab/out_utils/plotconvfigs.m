datadir = './data/conv/';
files = dir(datadir);


for i = 1:length(files)
    % init, skip ahead if needed
    f = files(i);
    if f.isdir
        continue;
    end
    name = f.name;
    
    % get conv, pe, th
    s = open([datadir,name]);
    
    cdiffs = s.conv.cdiffs;
    tlins = s.conv.tlins;
    terrs = s.conv.terrs;
    gnorms = s.conv.gnorms;
    idx = 0:(length(cdiffs)-1);
    
    
    figure;
    semilogy(idx,cdiffs,idx,tlins,idx,gnorms,idx,terrs,'LineWidth',2);
    hold on;
    legend('Class movement','Linear prob err','Rel grad norm','Testing err');
    xlabel('Newton iteration');
    ylabel('Newton metric score');
    title([name,' Newton metrics']);
    hold off
    
end