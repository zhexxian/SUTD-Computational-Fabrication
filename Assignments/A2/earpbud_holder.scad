module earbud_holder(height=5, earbud_radius=5, text_content = "ZX"){

    $fn=30;

    module base_top(){
        hull() {
            translate([12,0,0]) cylinder(h=height, r=earbud_radius);
            cylinder(h=height, r=earbud_radius);
        }
    }

    module base_middle(){
        minkowski() {
          cube([22,2,height]);
          sphere(2);
        }
    }

    module base_bottom(){
        difference(){
            hull () {
                cylinder(r=2,h=height+3);
                translate([10,0,0]) cylinder(r=2,h=height+3);
                translate([5,22,0]) cylinder(r=5,h=height+3);
            }
            {
                translate([-5, 22, 4]){
                    rotate(a=90,v=[0,1,0]){
                        cylinder(r=2.5,h=20);
                    }
                    
                }
                
                translate([-5, 20.5, 0]){
                    cube([20,3,height]);      
                }       
                
            }
        }
    }

    module top_divider(){
        hull () {
            cylinder(r=2,h=4);
            translate([5,0,0]) cylinder(r=2,h=4);
            translate([2.5,17,0]) cylinder(r=2,h=4);
        }
    }

    module earbud_container(){
        cylinder(h=10, r=earbud_radius);
    }

        
    module fingertip_container(){
        sphere(r=4);
    }

    module basic_shape(){
        translate([7,-10,2]){
            minkowski() {
               base_top();
               sphere(2);
            }   
        }

        translate([2,8,2]){
            base_middle();
        } 
        
        translate([8,-5,0]){
            minkowski() {
              base_bottom();
              sphere(1);
            }
        }

        translate([10.5,-8,-3]){
            minkowski() {
              top_divider();
              sphere(1);
            }
        }
    }


    module subtracted_shape(){

        for (x=[7, 19]) {
            translate([x,-10,-5]){
                earbud_container();
            }
        }
        
        for (x=[5, 21]) {
            translate([x,1,0]){
                fingertip_container();
            }
        }

    }

    difference(){
        difference(){

            basic_shape();
            subtracted_shape();
        }
        translate([-15,-25,height+2]){
            cube([50,50,5]);
        }
    }

    
    font = "Liberation Sans";

    translate ([11,21,-2]) {
        linear_extrude(height = 4) {
            rotate(a=180,v=[1,0,0]){
                text(text_content, font = font, size = 2);
            }
        }
    }
}
earbud_holder();
//earbud_holder(height=10, earbud_radius=4, text_content = "CF");