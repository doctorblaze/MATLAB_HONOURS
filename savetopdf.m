function savetopdf(filename,figurenumber)
    h = figure(figurenumber);
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,filename,'-dpdf','-r0')
end