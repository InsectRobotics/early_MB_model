This is a readme that briefly explains the use of AntData.mat for ant navigation studies. - Created by Mike Mangan 26/8/13

AntData is a struct that contains lots of informtation from the Seville, 2009 dataset.  For the purposes of ant navigation modelling.  The key data required are the 1cm spaced inward routes of ants, and also the headings along those routes.

These can be found in:

AntData.(AntX).InwardRouteData.(RouteY).One_cm_control_points

and 

AntData.(AntX).InwardRouteData.(RouteY).One_cm_control_points_headings

where X and Y must be defined by the user.  The control points are saved as [x,y] and the control point headings in degrees, where the positive X axis is 0 degrees, positive y axis is 90 degrees etc.

For early navigation studies I suggest we test with 3 ants only (Ant1, Ant2, Ant3), and with first four inward routes of each (Route1, Route2, Route3, Route4).

 
