function velbins_oct;
more off; warning off;
close all; clear all;

nidistzall=importdata('ni_distzall_527500000.dat');
nidistrall=importdata('ni_distrall_527470000.dat');

zlength = size(nidistzall)(1,1)/600; rheight = size(nidistrall)(1,1)/600;
velbinsr = zeros(1,600); velbinsz = zeros(1,600);
zbin_matrix = zeros(zlength,600); rbin_matrix = zeros(rheight,600);
tmpr = zeros(1,rheight); tmpz = zeros(1,zlength);

 for vbin=1:600;
 	for z=1:zlength;
	
		normvz = sum ( nidistzall ( (z-1)*600+1:1:z*600, 3 ) ) / 600;
		% printf('>> normvz is %d\n', normvz);
		
		if ( normvz == 0) 
			zbin_matrix ( z, vbin ) = 0; 
		else
 			tmpz ( z ) = sum ( sum ( nidistzall ( (z-1)*600+vbin:600:end , 4:6 ) ...
 									 ./ normvz ) );
 			velbinsz ( vbin ) += tmpz ( z );
 			zbin_matrix ( z, vbin ) = tmpr ( z ); 
		end

 	end

	for r=1:rheight;
	
		normvr = sum ( nidistrall ( (r-1)*600+1:1:r*600, 3 ) ) / 600;
		% printf('>> normvr is %d\n', normvr);
		
		if ( normvr == 0)
			velbinsr ( vbin ) += 0;
			rbin_matrix ( r, vbin ) = 0; 
		else
 			tmpr ( r ) = sum ( sum ( nidistrall ( (r-1)*600+vbin:600:end , 4:6 ) ...
 									 ./ normvr ) );
 			velbinsr ( vbin ) += tmpr ( r );
 			rbin_matrix ( r, vbin ) = tmpr ( r ); 
		end

 	end

end

keyboard;
save -text 'tmp.dat' *
clear all;

end
