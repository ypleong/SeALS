
%> directory file
file = '.';
% Call the documentacion with standard settings (to run the graph GraphViz
% has to be installed)
m2html('mfiles', file , 'htmldir', 'doc', 'recursive','on', 'global','on', ...
      'template','frame', 'index','menu', 'graph','on','save','on');
