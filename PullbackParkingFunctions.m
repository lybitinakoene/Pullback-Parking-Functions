~}%% Main Function
function [] = PullbackParkingFunctions()

	%% Introduction
	clc
    dbstop if error
    dbstop if naninf


	% Iterates all functions from 1 to this constant
	HIGH_N = 6;


	% Determines which equations/algorithms to run
	RECURSIVE_FORMULA = 1;
	BRUTE_FORCE_ALGORITHM = 1;
	COUNTING_THROUGH_PERMUTATIONS = 1;



	%% Formulaic iteration
	if RECURSIVE_FORMULA

		% Allocates space for the variables to present at the end of the
        % function
		recursiveTotals = zeros(HIGH_N, HIGH_N, HIGH_N, HIGH_N);
		recursiveCaseTotals = zeros(HIGH_N, HIGH_N, 3, HIGH_N, HIGH_N);

		% Sets the base case of m=1, n=1 for all k and l
		for i = 1:HIGH_N
            for j = 1:HIGH_N
			    recursiveTotals(1, 1, i, j) = 1;
			    recursiveCaseTotals(1, 1, 1, i, j) = 1;
            end
		end

		% Iterates over each PF_m, n(k, l)
		for n = 2:HIGH_N
			for m = 1:n
				for k = 0:HIGH_N-1
                    for l = 0:HIGH_N-1

					    % Iterates over each possible spot for the mth
                        % car to end up in
					    for i = 1:n

						    % Calculates case 2a (this is only
                            % calculated for x=0, so it's outside of the
                            % x iteration)
						    case2ATotals = case2AFunction(...
                                recursiveTotals, n, m, k, l, i);
						    recursiveTotals(m, n, k+1, l+1) = ...
                                recursiveTotals(m, n, k+1, l+1) + ...
                                case2ATotals;
						    recursiveCaseTotals(m, n, 2, k+1, l+1) = ...
                                recursiveCaseTotals(...
                                m, n, 2, k+1, l+1) + case2ATotals;

						    % Calculates case 3a (this is only
                            % calculated for x=0, so it's outside of the
                            % x iteration)
						    case3ATotals = case3AFunction(...
                                recursiveTotals, n, m, k, l, i);
						    recursiveTotals(m, n, k+1, l+1) = ...
                                recursiveTotals(m, n, k+1, l+1) + ...
                                case3ATotals;
						    recursiveCaseTotals(m, n, 3, k+1, l+1) = ...
                                recursiveCaseTotals(...
                                m, n, 3, k+1, l+1) + case3ATotals;

						    % Iterates over all posible choices for
                            % number of cars on the leftmost side of the
                            % street
						    for x = 0:m-1

							    % Calculates case 1
							    case1Totals = case1Function(...
                                    recursiveTotals, n, m, k, l, i, x);
							    recursiveTotals(m, n, k+1, l+1) = ...
                                    recursiveTotals(m, n, k+1, l+1) +...
                                    case1Totals;
							    recursiveCaseTotals(...
                                    m, n, 1, k+1, l+1) = ...
                                    recursiveCaseTotals(m, n, 1, k+1, ...
                                    l+1) + case1Totals;

							    % Calculates case 2
							    case2BTotals = case2BFunction(...
                                    recursiveTotals, n, m, k, l, i, x);
							    recursiveTotals(m, n, k+1, l+1) = ...
                                    recursiveTotals(m, n, k+1, l+1) +...
                                    case2BTotals;
							    recursiveCaseTotals(...
                                    m, n, 2, k+1, l+1) = ...
                                    recursiveCaseTotals(m, n, 2, k+1, ...
                                    l+1) + case2BTotals;

							    % Calculates case 3
							    case3BTotals = case3BFunction(...
                                    recursiveTotals, n, m, k, l, i, x);
							    recursiveTotals(m, n, k+1, l+1) = ...
                                    recursiveTotals(m, n, k+1, l+1) +...
                                    case3BTotals;
							    recursiveCaseTotals(...
                                    m, n, 3, k+1, l+1) = ...
                                    recursiveCaseTotals(m, n, 3, k+1, ...
                                    l+1) + case3BTotals;

						    end
					    end

                    end
				end
			end
		end

	end



	%% Brute force iteration
	if BRUTE_FORCE_ALGORITHM

		% Allocates space for the variables to present at the end of the
        % function
		bruteForceTotals = zeros(HIGH_N, HIGH_N, HIGH_N, HIGH_N);
		bruteForceCaseTotals = zeros(HIGH_N, HIGH_N, 3, HIGH_N, HIGH_N);

		% Iterates over all m and n
		for n = 1:HIGH_N
			for m = 1:n

				% Defines initial preference vector of all ones except
                % the final spot
				preferenceVector = ones(m+1, 1);
				preferenceVector(m+1) = 0;

				% Iterates over all possible preference vectors
				while 1

					% Iterates over all k and l
					for k = 0:HIGH_N - 1
                        for l = 0:HIGH_N - 1

						    % Calculates which case, if any, the
                            % preference vector represents.  All
                            % variables will be 0 if the preference
                            % vector doesn't park all cars
						    [case1, case2, case3] = bruteForcePark(...
                                preferenceVector, n , m , k, l);

						    % Adds one to the correct case
                            if case1
							    bruteForceTotals(m, n, k+1, l+1) = ...
                                    bruteForceTotals(...
                                    m, n, k+1, l+1) + 1;
							    bruteForceCaseTotals(...
                                    m, n, 1, k+1, l+1) = ...
                                    bruteForceCaseTotals(m, n, 1, ...
                                    k+1, l+1) + 1;
						    elseif case2
							    bruteForceTotals(m, n, k+1, l+1) = ...
                                    bruteForceTotals(...
                                    m, n, k+1, l+1) + 1;
							    bruteForceCaseTotals(...
                                    m, n, 2, k+1, l+1) = ...
                                    bruteForceCaseTotals(m, n, 2, ...
                                    k+1, l+1) + 1;
						    elseif case3
							    bruteForceTotals(m, n, k+1, l+1) = ...
                                    bruteForceTotals(...
                                    m, n, k+1, l+1) + 1;
							    bruteForceCaseTotals(...
                                    m, n, 3, k+1, l+1) = ...
                                    bruteForceCaseTotals(m, n, 3, ...
                                    k+1, l+1) + 1;
                            end

                        end
					end

					% Chooses the next preference vector
                    for i = 1:m+1
						preferenceVector(i) = preferenceVector(i) + 1;
						if preferenceVector(i) ~~= n+1
							break
						else
							preferenceVector(i) = 1;
						end
                    end

					% If all preference vectors have been tested, break out
                    % of the loop
					if preferenceVector(m+1) == 1
						break;
					end

				end
			end
		end

	end



	%% Counting through permutations formula
	if COUNTING_THROUGH_PERMUTATIONS

		% Allocates space for the variables to present at the end of the
        % function
		countingPermutationsTotals = zeros(HIGH_N, HIGH_N, HIGH_N, HIGH_N);

		% Iterates over all m, n, and k
		for n = 1:HIGH_N
			for m = 1:n

				% Calculates all permutations of m and n recursively
				[permList, permListLength] = getPermutations(1, ...
                    zeros(n, 1), m, n);

				for k = 0:HIGH_N-1
                    for l = 0:HIGH_N-1

					    % Sums over all permutations
					    for piIndex = 1:permListLength

						    % Multiplies all possibilities of each spot
                            % together
						    permProduct = 1;
						    for i = 1:n

							    possiblePrefsNumber = functionB(i, ...
                                    permList(:, piIndex), m, n, k, l) + ...
                                    functionF(i, permList(:, piIndex), ...
                                    m, n, k, l) + 1;
							    permProduct = permProduct * ...
                                    possiblePrefsNumber;

						    end
						    countingPermutationsTotals(...
                                m, n, k+1, l+1) = ...
                                countingPermutationsTotals(m, n, ...
                                k+1, l+1) + permProduct;

					    end

                    end
				end
			end
		end

	end



    %% Sa
    save("KLMNOP.mat")



	%% Pause to look at found data
    if RECURSIVE_FORMULA && BRUTE_FORCE_ALGORITHM && ...
            COUNTING_THROUGH_PERMUTATIONS
        for n = 1:HIGH_N
            for m = 1:HIGH_N
                for k = 1:HIGH_N
                    for l = 1:HIGH_N

                        if bruteForceTotals(m, n, k, l) ~~= recursiveTotals(...
                                m, n, k, l)
                            NaN;
                        end

                        if bruteForceTotals(m, n, k, l) ~~= ...
                                countingPermutationsTotals(m, n, k, l)
                            NaN;
                        end

                        for i = 1:3
                            if bruteForceCaseTotals(m, n, i, k, l) ~~= ...
                                    recursiveCaseTotals(m, n, i, k, l)
                                NaN;
                            end
                        end

                    end
                end
            end
        end
    end
	NaN;

end


%% getPermutations
% Calculates all permutations of m cars and n spots recursively
%
% mIndex is the car we are parking during this iteration
% nList is the list of all spots which are still open
% m is total number of cars
% n is total number of spots
%
% permList is the full list of possible permutations
% listLength is the length of permList
%
% % % % %
function [permList, listLength] = getPermutations(mIndex, nList, m, n)

    % Defines the length of the current iteration of lists, i.e. all
    % possibilities for parking the remaining cars in the remaining spots
    listLength = 0;

    % Iterates over all of the spots
    for i = 1:n

        % If there are no cars left, the currently parked cars is the only
        % possibility, so return this permutation
        if mIndex == m+1
            permList(:, 1) = nList;
            listLength = 1;
            break
        end

        % Insures that the spot is open
        if nList(i) == 0

            % Park car mIndex in spot i
            nList(i) = mIndex;

            % Recurse over all possibilities for this setup
            [nextPermList, nextListLength] = ...
                getPermutations(mIndex+1, nList, m, n);

            % Adds these new permutations to the list and increases
            % length accordingly
            permList(:, listLength+1:listLength+nextListLength) = ...
                nextPermList;
            listLength = listLength + nextListLength;

            % Resets the spot to 0 for the next i iteration
            nList(i) = 0;

        end

    end

end


%% functionB
% Calculates B (number of possible preferences car i could have had to the
% right of its chosen spot) for a given spot and permutation.  B is defined
% as B = min(R, k)
%
% i is the current car we're looking at
% pi is the current permutation
% m is total number of cars
% n is total number of spots
% k is the amount the car is permitted to back up
% l is the amount the car is permitted to drive forward
%
% bTotal is B
%
% % % % %
function [bTotal] = functionB(i, pi, m, n, k, l)

	bTotal = min(functionR(i, pi, m, n, k, l), k);

end


%% functionF
% Calculates F (number of possible preferences car i could have had to the
% left of its chosen spot) for a given spot and permutation.  F is defined
% as F = min(i-1, l) when L = i-1, and defined as F = max(min(L-k, l), 0)
% otherwise
%
% i is the current car we're looking at
% pi is the current permutation
% m is total number of cars
% n is total number of spots
% k is the amount the car is permitted to back up
% l is the amount the car is permitted to drive forward
%
% fTotal is F
%
% % % % %
function [fTotal] = functionF(i, pi, m, n, k, l)

    variableL = functionL(i, pi, m, n, k, l);
    if variableL == i - 1
        fTotal = min(i-1, l);
    else
	    fTotal = max(min(variableL - k, l), 0);
    end

end


%% functionFContainedVersion
% Calculates F (number of possible preferences car i could have had to the
% left of its chosen spot) for a given spot and permutation for the
% contained counting through permutations function.  F is defined as
% F = max(min(L-k, l), 0)
%
% i is the current car we're looking at
% pi is the current permutation
% m is total number of cars
% n is total number of spots
% k is the amount the car is permitted to back up
% l is the amount the car is permitted to drive forward
%
% fTotal is F
%
% % % % %
function [fTotal] = functionFContainedVersion(i, pi, m, n, k, l)

    fTotal = max(min(functionL(i, pi, m, n, k, l) - k, l), 0);

end


%% functionR
% Calculates the number of spots to the right of car i which are taken
%
% i is the current car we're looking at
% pi is the current permutation
% n is total number of spots
%
% rTotal is R
%
% % % % %
function [rTotal] = functionR(i, pi, ~~, n, ~~, ~~)

    % Iterates over all spots right of spot i
	rTotal = 0;
	for spot = i+1:n

        % If the spot is filled with a car that would have parked before
        % the car in spot i, add one to total and continue, else break out
        % of loop
        if pi(spot) < pi(i) && pi(spot) ~~= 0
			rTotal = rTotal + 1;
		else
			break
        end

	end

end


%% functionL
% Calculates the number of spots to the left of car i which are taken
%
% i is the current car we're looking at
% pi is the current permutation
%
% lTotal is L
%
% % % % %
function [lTotal] = functionL(i, pi, ~~, ~~, ~~, ~~)

    % Iterates over all spots left of spot i
	lTotal = 0;
	for spot = i-1:-1:1

        % If the spot is filled with a car that would have parked before
        % the car in spot i, add one to total and continue, else break out
        % of loop
        if pi(spot) < pi(i) && pi(spot) ~~= 0
			lTotal = lTotal + 1;
		else
			break
        end

	end

end


%% containedParking
% Calculates the number of possible preference lists for a contained
% k, l, m, n parking function. When we are not considering l, this is a
% very simple equation, but in this case, k, l, m, n means we need to run
% this process using counting through permutations
%
% m is total number of cars
% n is total number of spots
% k is the amount the car is permitted to back up
% l is the amount the car is permitted to drive forward
%
% total is |B~~_{m, n}|
%
% % % % %
function [total] = containedParking(m, n, k, l)

	% Calculates all permutations of m and n recursively
    [permList, permListLength] = getPermutations(1, zeros(n, 1), m, n);
    total = 0;

    % Sums over all permutations
    for piIndex = 1:permListLength

	    % Multiplies all possibilities of each spot together
	    permProduct = 1;
	    for i = 1:n

            possiblePrefsNumber = functionB(i, permList(:, piIndex), ...
                m, n, k, l) + functionFContainedVersion(...
                i, permList(:, piIndex), m, n, k, l) + 1;
		    permProduct = permProduct * possiblePrefsNumber;

	    end
	    total = total + permProduct;

    end

end


%% Case functions
% For each case below, the calculations should match the formula given in
% the "Pullback Parking Functions" paper (Elder et al.) precisely, but if
% you find a discrepancy, please report it to one of the authors.  The
% insides of the functions will not be commented because most of it is very
% similar to the rest, and follows the logic described in the paper.  There
% are four different logical parts used in the cases over and over, and I
% will describe these below.
%
% Choose functions (nchoosek) are used throughout the cases to choose a set
% of cars from the remaining cars.  It uses if-else logic to ensure that
% the choose function only allows the correct parameters, and otherwise
% will say that this case is not possible, and return 0.
%
% Callbacks to the previously calculated parking function totals
% (recursiveTotals) are also used throughout.  It will return the correct
% number if the indices are valid within the recursiveTotals matrix.
% Otherwise, if 0 <= mVal <= nVal, then it will return 1, otherwise, 0.
%
% Calls to the Contained Parking Function totals (containedParking) are
% also used.  This function works using the same logic as counting through
% permutations, and so works for 0 < mVal <= nVal.  We need to additionally
% consider the case mVal = 0 and mVal <= nVal, in which case we should
% return 1.  Otherwise, we will return 0.
%
% Occasionally, the total is multiplied by a min function.  The reasons for
% this vary by use case, and do not need any if-else statements as the
% definition for the min function allows for any two inputs.


%% case1Function
% Calculates totals for case 1 given the inputs as parameters
%
% recursiveTotals is the matrix containing the totals for previously
% calculated recursive values
% n is total number of spots
% m is total number of cars
% k is the amount the car is permitted to back up
% l is the amount the car is permitted to drive forward
% i is the spot which car m ends up parking in (in this case, also its
% preference)
% x is the number of cars which will park to the left of car m
%
% caseTotals is the calculated total for this case
%
% % % % %
function [caseTotals] = case1Function(recursiveTotals, n, m, k, l, i, x)

	if m-1 >= x && x >= 0
		caseTotals = nchoosek(m-1, x);
	else
		caseTotals = 0;
	end

	if i-1 > 0 && x > 0
		caseTotals = caseTotals * recursiveTotals(x, i-1, k+1, l+1);
	elseif i-1 < x || x < 0
		caseTotals = 0;
	end

    if n-i >= m-1-x && m-1-x > 0
		caseTotals = caseTotals * containedParking(m-1-x, n-i, k, l);
    elseif n-i < m-1-x || m-1-x < 0
		caseTotals = 0;
    end

end


%% case2BFunction
% Calculates totals for case 2b given the inputs as parameters
%
% recursiveTotals is the matrix containing the totals for previously
% calculated recursive values
% n is total number of spots
% m is total number of cars
% k is the amount the car is permitted to back up
% l is the amount the car is permitted to drive forward
% i is the spot which car m ends up parking in
% x is the number of cars which will park to the left of car m
%
% caseTotals is the calculated total for this case
%
% % % % %
function [caseTotals] = case2BFunction(recursiveTotals, n, m, k, l, i, x)

	if m-1 >= x && x >= 0
		caseTotals = nchoosek(m-1, x);
	else
		caseTotals = 0;
	end

	if i-1 > 0 && x > 0
		caseTotals = caseTotals * recursiveTotals(x, i-1, k+1, l+1);
	elseif i-1 < x || x < 0
		caseTotals = 0;
	end

	helperTotals = 0;
	for R = 1:n-i-1
        helperTotals = helperTotals + case2FunctionHelper(...
            recursiveTotals, n, m, k, l, i, x, R);
	end
	caseTotals = caseTotals * helperTotals;

end


%% case2FunctionHelper
% Calculates totals for case 2b's iteration over R given the inputs as
% parameters
%
% n is total number of spots
% m is total number of cars
% k is the amount the car is permitted to back up
% l is the amount the car is permitted to drive forward
% i is the spot which car m ends up parking in
% x is the number of cars which will park to the left of car m
% R is the number of cars which park in the "populated region", in this
% case, the densely populated region to the right of car m
%
% caseTotals is the calculated total for this case
%
% % % % %
function [caseTotals] = case2FunctionHelper(~~, n, m, k, l, i, x, R)

	if m-1-x >= R && R >= 0
		caseTotals = nchoosek(m-1-x, R);
	else
		caseTotals = 0;
	end

	if R > 0
		caseTotals = caseTotals * containedParking(R, R, k, l);
    elseif R < 0
		caseTotals = 0;
	end

	if n-R-i-1 >= m-1-x-R && m-1-x-R > 0
		caseTotals = caseTotals * containedParking(m-1-x-R, n-R-i-1, k, l);
    elseif n-R-i-1 < m-1-x-R || m-1-x-R < 0
		caseTotals = 0;
	end

	caseTotals = caseTotals * min(R, k);

end


%% case2AFunction
% Calculates totals for case 2a given the inputs as parameters
%
% recursiveTotals is the matrix containing the totals for previously
% calculated recursive values
% n is total number of spots
% m is total number of cars
% k is the amount the car is permitted to back up
% l is the amount the car is permitted to drive forward
% i is the spot which car m ends up parking in
%
% caseTotals is the calculated total for this case
%
% % % % %
function [caseTotals] = case2AFunction(recursiveTotals, n, m, k, l, i)

	if m-1 >= n-i && n-i >= 0
		caseTotals = nchoosek(m-1, n-i);
	else
		caseTotals = 0;
	end

    if i-1 > 0 && m-1-n+i > 0
	    caseTotals = caseTotals * recursiveTotals(m-1-n+i, i-1, k+1, l+1);
	elseif m-1-n+i > i-1 || i-1 < 0
		caseTotals = 0;
    end

    if n-i > 0
	    caseTotals = caseTotals * containedParking(n-i, n-i, k, l);
    elseif n-1 < 0
        caseTotals = 0;
    end

	caseTotals = caseTotals * min(k, n-i);

end


%% case3BFunction
% Calculates totals for case 3b given the inputs as parameters
%
% recursiveTotals is the matrix containing the totals for previously
% calculated recursive values
% n is total number of spots
% m is total number of cars
% k is the amount the car is permitted to back up
% l is the amount the car is permitted to drive forward
% i is the spot which car m ends up parking in
% x is the number of cars which will park to the left of car m
%
% caseTotals is the calculated total for this case
%
% % % % %
function [caseTotals] = case3BFunction(recursiveTotals, n, m, k, l, i, x)

	if m-1 >= x && x >= 0
		caseTotals = nchoosek(m-1, x);
	else
		caseTotals = 0;
	end

	helperTotals = 0;
    for R = k+1:i-2
        helperTotals = helperTotals + case3FunctionHelper(...
            recursiveTotals, n, m, k, l, i, x, R);
    end
	caseTotals = caseTotals * helperTotals;

end


%% case3FunctionHelper
% Calculates totals for case 3b's iteration over R given the inputs as
% parameters
%
% recursiveTotals is the matrix containing the totals for previously
% calculated recursive values
% n is total number of spots
% m is total number of cars
% k is the amount the car is permitted to back up
% l is the amount the car is permitted to drive forward
% i is the spot which car m ends up parking in
% x is the number of cars which will park to the left of car m
% R is the number of cars which park in the "populated region", in this
% case, the densely populated region to the right of car m
%
% caseTotals is the calculated total for this case
%
% % % % %
function [caseTotals] = case3FunctionHelper(recursiveTotals, ...
    n, m, k, l, i, x, R)

    if i-R-2 > 0 && x > 0
		caseTotals = recursiveTotals(x, i-R-2, k+1, l+1);
    elseif i-R-2 < x || x < 0
		caseTotals = 0;
    else
        caseTotals = 1;
    end

	if m - 1 - x >= R && R >= 0
		caseTotals = caseTotals * nchoosek(m-1-x, R);
	else
		caseTotals = 0;
	end

	if R > 0
		caseTotals = caseTotals * containedParking(R, R, k, l);
    elseif R < 0
		caseTotals = 0;
	end

	if n-i >= m-1-x-R && m-1-x-R > 0
		caseTotals = caseTotals * containedParking(m-1-x-R, n-i, k, l);
    elseif n-i < m-1-x-R || m-1-x-R < 0
		caseTotals = 0;
	end

	caseTotals = caseTotals * min(R-k, l);

end


%% case3AFunction
% Calculates totals for case 3a given the inputs as parameters
%
% recursiveTotals is the matrix containing the totals for previously
% calculated recursive values
% n is total number of spots
% m is total number of cars
% k is the amount the car is permitted to back up
% l is the amount the car is permitted to drive forward
% i is the spot which car m ends up parking in
%
% caseTotals is the calculated total for this case
%
% % % % %
function [caseTotals] = case3AFunction(recursiveTotals, n, m, k, l, i)

	if m-1 >= i-1 && i-1 >= 0
		caseTotals = nchoosek(m-1, i-1);
	else
		caseTotals = 0;
	end

    if i-1 > 0
	    caseTotals = caseTotals * recursiveTotals(i-1, i-1, k+1, l+1);
	elseif i-1 < 0
		caseTotals = 0;
    end

    if k == 0 && l == 1 && m == 5 && n == 5 && i == 2
        NaN;
    end

    if n-i >= m-i && m-i > 0
	    caseTotals = caseTotals * containedParking(m-i, n-i, k, l);
    elseif n-i < m-i || m-i < 0
        caseTotals = 0;
    end

	caseTotals = caseTotals * min(i-1, l);

end


%% bruteForcePark
% Finds out if a given preference vector would be a valid parking function
%
% preferenceVector is the preference vector to check
% n is total number of spots
% m is total number of cars
% k is the amount the car is permitted to back up
% l is the amount the car is permitted to drive forward
%
% case1 is a boolean for the case where car m gets its preference
% case2 is a boolean for the case where car m backs into its spot
% case3 is a boolean for the case where car m moves forward into its spot
%
% % % % %
function [case1, case2, case3] = bruteForcePark(preferenceVector, ...
    n, m, k, l)

    % Sets the cases to zero.  All cases should return 0 if the preference
    % vector is not a valid parking function, and otherwise only one will
    % return 1.
	case1 = 0;
	case2 = 0;
	case3 = 0;

    % The variable for car m's spot
	finalSpot = 0;

    % Set the street to be empty
	parkingSpots = zeros(n, 1);

    % Iterates over all cars
    for i = 1:m

        carHasParked = 0;

        % If the car's preference is open, park there, otherwise try
        % backing up and then pulling forward
        if parkingSpots(preferenceVector(i)) == 0

            % Park the car in its preferred spot
            parkingSpots(preferenceVector(i)) = 1;
            carHasParked = 1;

            % If the car is the last car, set the finalSpot variable
			if i == m
				finalSpot = preferenceVector(i);
			end

        else

            % Iterates over all spots which are up to k behind the car's
            % preference
            index = 1;
            currentSpot = preferenceVector(i) - index;
            while index <= k && currentSpot > 0

                % If the current spot is open
                if parkingSpots(currentSpot) == 0

                    % Park the car
                    parkingSpots(currentSpot) = 1;
                    carHasParked = 1;

                    % If the car is the last car, set the finalSpot
                    % variable
                    if i == m
						finalSpot = currentSpot;
                    end

                    % Break out of the loop so the car only parks once
                    break;

                end

                % Sets the next spot to check
                index = index + 1;
                currentSpot = preferenceVector(i) - index;

            end

            % If the car has still not parked after backing up
            if ~~carHasParked

                % Iterates over all spots in front of car i's preference
                % up to i + l
                index = 1;
                currentSpot = preferenceVector(i) + index;
                while currentSpot < n+1 && index <= l

                    % If the spot is free
                    if parkingSpots(currentSpot) == 0

                        % Park the car
                        parkingSpots(currentSpot) = 1;
                        carHasParked = 1;

                        % If the car is the last car, set the finalSpot
                        % variable
                        if i == m
							finalSpot = currentSpot;
                        end

                        % Break out of the loop so the car only parks once
                        break;

                    end

                    % Sets the next spot to check
                    index = index + 1;
                    currentSpot = preferenceVector(i) + index;

                end
            end

        end

        % If the car wasn't able to park, break out of the loop
        if ~~carHasParked
            break
        end

    end

    % If all cars parked
	if finalSpot ~~= 0

        % Set the correct case
        if finalSpot == preferenceVector(m)
            case1 = 1;
        elseif finalSpot < preferenceVector(m)
            case2 = 1;
        elseif finalSpot > preferenceVector(m)
            case3 = 1;
        end

	end

end


