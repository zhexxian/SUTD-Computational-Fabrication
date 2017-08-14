module flexureSolid2D(width, height, wCut, hStalk) {

	//cube from which flexure will be cut
	hEnd = 0.5*(height - hStalk);
	wStalk = width - 2*wCut;

	difference() {
		difference() {
			square(size=[width, height], center=true);
			translate([-(wCut + wStalk)/2, 0, 0]) {
				square(size=[wCut, hStalk], center=true);
			}
		}

		translate([(wCut + wStalk)/2, 0, 0]) {
			square(size=[wCut, hStalk], center=true);
		}
	}	
}

module flexure(wCut, thickness, depth)
{
	width = 20;
    thickness = 20*thickness;
    depth = 20*depth;
    
	height = width;
	wCut = width*wCut;
	thickness = width*thickness;
	thickness2 = thickness;
	cavitySize = (0.7*height - 2*thickness)/4;
	hStalk2 = 2*cavitySize;

	echo("Thickness ", thickness2);
	echo("Stalk Height ", hStalk2);

	//deal with constraints
	linear_extrude(height = depth, center = true) {
		difference() {
			flexureSolid2D(width, height, wCut, hStalk2);
			flexureSolid2D(width-2*thickness, 0.7*height, wCut, hStalk2+2*thickness2);
		}
	}
}


flexure(0.25,0.05,0.4); 