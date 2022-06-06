function DoMyViz_node_rel(cortex, idx_atlas, var_node, filename, info)

    close all;
    N = length(cortex.Atlas(idx_atlas).Scouts);
    col = zeros(N,3);

    k1=figure(1)
    p=var_node; 
    val_max=max(p);
    val_min=min(p);

    for i=1:N

         nver = length(cortex.Atlas(idx_atlas).Scouts(i).Vertices);

         col(cortex.Atlas(idx_atlas).Scouts(i).Vertices,:) = repmat(cmap(p(i)/val_max,'summer'),nver,1);

    end
    colmap = cmap(64,'summer');
    set(gcf,'color','w')
    hp2 = patch('faces',cortex.Faces,'vertices',cortex.Vertices,'FaceVertexCData',col);

    set(hp2,'linestyle','none','FaceColor', 'interp','specularstrength',0);
    axis equal tight off
    alpha(1)
    light('position',[-1 0 0]);
    light('position',[1 0 0]);
    light('position',[0 0 1]);
    lighting gouraud
    view([-90 90]); 
    title(info, 'Interpreter', 'none');
    colormap(colmap);
    un=unique(p);
    ticks{1}=num2str(val_min);
    ticks{2}=num2str(val_max);
    c = colorbar('TickLabels',ticks, ...
                   'Ticks', [0,  1.0]);
    c.Title.String='Reliability';
%     saveas(k1,filename);
%     saveas(gcf,strcat(filename,'.pdf'));
    print(gcf, strcat(filename,'.png'), '-dpng','-r600');


end