function[Ind_keep] = Remove_irrelevant_data(Data)

% Remove irrelevant intersection points

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

% input:
% - Data: loaded raw data (result of ray tracing procedure)

% output: 
% - row indices of signal samples carrying relevant information

% Find non-zero amplitudes --> have to be kept
Ind_non_zero_amp = find(Data(:,4) ~= 0);

%--------------------------------------------------------------------------

% Indices: fivefold + zero amplitude
Ind_1 = find(Data(:,5) == 5);
Ind_five = intersect(Ind_non_zero_amp,Ind_1);
Ind_five_1 = Ind_five-1;
Ind_five_2 = Ind_five-2;
Ind_five_3 = Ind_five-3;
Ind_five_4 = Ind_five-4;

Ind_five = union(Ind_five,Ind_five_1);
Ind_five = union(Ind_five,Ind_five_2);
Ind_five = union(Ind_five,Ind_five_3);
Ind_five = union(Ind_five,Ind_five_4);

% Union
Ind_keep = union(Ind_non_zero_amp,Ind_five);

%--------------------------------------------------------------------------

% Indices: fourfold + zero amplitude
Ind_2 = find(Data(:,5) == 4);
Ind_four = intersect(Ind_non_zero_amp,Ind_2);
Ind_four_1 = Ind_four-1;
Ind_four_2 = Ind_four-2;
Ind_four_3 = Ind_four-3;

Ind_four = union(Ind_four,Ind_four_1);
Ind_four = union(Ind_four,Ind_four_2);
Ind_four = union(Ind_four,Ind_four_3);

% Union
Ind_keep = union(Ind_keep,Ind_four);

%--------------------------------------------------------------------------

% Indices: triple + zero amplitude
Ind_3 = find(Data(:,5) == 3);
Ind_trip = intersect(Ind_non_zero_amp,Ind_3);
Ind_trip_1 = Ind_trip-1;
Ind_trip_2 = Ind_trip-2;

Ind_trip = union(Ind_trip,Ind_trip_1);
Ind_trip = union(Ind_trip,Ind_trip_2);

% Union
Ind_keep = union(Ind_keep,Ind_trip);
%--------------------------------------------------------------------------

% Indices: double + zero amplitude
Ind_4 = find(Data(:,5) == 2);
Ind_doub = intersect(Ind_non_zero_amp,Ind_4);
Ind_doub_1 = Ind_doub-1;

Ind_doub = union(Ind_doub,Ind_doub_1);

% Union
Ind_keep = union(Ind_keep,Ind_doub);

