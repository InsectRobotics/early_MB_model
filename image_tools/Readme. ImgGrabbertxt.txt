Readme for ImgGrabber.m function.  Readme written by Mike Mangan 30/8/13.

The ImgGrabber function grabs either "robot-eye" or "ant-eye" images from the world500_gray.

Current settings for robot eye are hfov=360 and resolution=1.
Current settings for ant eye are hfov=296 and resolution=4.

Note: the function has yet to be tested in other settings so stick to these settings for now.

The images produced by the robot eye are 75*360 in size, where the vertical axis encodes elevations from -14.5 degrees to 59.5 degrees in 1 degree steps.  Horizontal pixals encode sampling from -179.5 to 179.5 in 1 degree intervals.  This has been confirmed using a toy world with cylinders of known height, width and distance.

The images produced by the ant eye are 19*74 in size, where the vertical axis encodes elevations from -14.5 degrees to 57.5 degrees in 4 degree steps.  Horizontal pixals encode sampling from -146 to 146 in 4 degree intervals.  The images are firstly blurred at 2degrees using the blockproc procedure to be like the real ant acceptance and interommnitidial angles.
