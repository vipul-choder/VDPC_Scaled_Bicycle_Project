function xdot = scaled_bicycle_model(t,x)

global Mphiphi Mdelphi Mdeldel Mphidel

global Kophiphi Kophidel Kodelphi Kodeldel K2phiphi K2phidel K2delphi K2deldel

global Cophiphi Cophidel Codelphi Codeldel 

global v gravacc

roll_angle = x(1);
roll_rate = x(2);
steer_angle = x(3);
steer_rate = x(4);

%Addition for Step 2

xdot(1,:) = roll_rate;

xdot(2,:) = (1/(Mphidel*Mdelphi - Mdeldel*Mphiphi))*(roll_rate*v*(Cophiphi*Mdelphi*Mdeldel - Codelphi*Mphidel*Mdelphi)...
                            + steer_rate*v*(Mdeldel*Cophidel*Mdelphi - Codeldel*Mphidel*Mdelphi)...
                            + roll_angle*(gravacc*Kophiphi*Mdelphi*Mdeldel + v^2*Mdeldel*K2phiphi*Mdelphi - gravacc*Mdeldel*Kodelphi*Mphiphi - v^2*K2delphi*Mphiphi*Mdeldel -gravacc*Kodelphi*Mphidel*Mdelphi + gravacc*Kodelphi*Mdeldel*Mphiphi - v^2*K2delphi*Mphidel*Mdelphi + v^2*K2delphi*Mdeldel*Mphiphi)...
                            + steer_angle*(gravacc*Mdeldel*Kodelphi*Mdelphi + v^2*Mdeldel*K2phidel*Mdelphi - gravacc*Mdeldel*Kodeldel*Mphiphi - v^2*Mdeldel*K2deldel*Mphiphi - gravacc*Kodeldel*Mdelphi*Mphidel + gravacc*Kodeldel*Mdeldel*Mphiphi- v^2*K2deldel*Mphidel*Mdelphi + v^2*K2deldel*Mdeldel*Mphiphi));% rolling acceleration

xdot(3,:) = steer_rate;

xdot(4,:) = (1/(Mphidel*Mdelphi - Mdeldel*Mphiphi))*(roll_rate*v*(-Cophiphi*Mdelphi + Codelphi*Mphiphi)...
                                                    + steer_rate*v*(Codeldel*Mphiphi - Cophidel*Mdelphi)...
                                                    - steer_angle*(gravacc*Kophidel*Mdelphi + v^2*K2phidel*Mdelphi - gravacc*Kodeldel*Mphiphi - v^2*K2deldel*Mphiphi)...
                                                    - roll_angle*(gravacc*Kophiphi*Mdelphi + v^2*K2phiphi*Mdelphi - gravacc*Kodelphi*Mphiphi - v^2*K2delphi*Mphiphi));
%----------------------------------------------------------------
