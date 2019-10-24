#include "PatternMatcher.h"
#include "utils.h"

using namespace std;

vector<int> prepareLPS(vector<char> pattern);

vector<int> PatternMatcher::naive()
{
	vector<int> matches;

	size_t N = text.size();
	size_t M = pattern.size();

	if (N > 0 && M > 0 && M <= N)
	{
		Timer t;
		t.start();

		for (int i = 0; i < N - M + 1; i++) {
			int j = 0;

			while (j < M && text[i + j] == pattern[j]) j++;

			if (j == M)
			{
				matches.push_back(i);
			}
		}
		t.stop();

		cout << "[NAIVE MATCH]: Text length: " << text.size() <<
			", Pattern length: " << pattern.size() <<
			" -> Time: " << t.elapsed() << endl;
	}
	else
	{
		cout << "[NAIVE MATCH]: Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::boyer_moore()
{
	vector<int> matches;

	size_t N = text.size();
	size_t M = pattern.size();

	if (N > 0 && M > 0 && M <= N)
	{
		// TODO
	}
	else
	{
		cout << "[BOYER-MOORE MATCH] Text and pattern did not load properly...";
	}

	return matches;
}

vector<int> PatternMatcher::kmp()
{
	vector<int> matches;

	size_t N = text.size();
	size_t M = pattern.size();

	if (N > 0 && M > 0 && M <= N)
	{
		Timer t;
		t.start();

		vector<int> lps(M, 0);

		lps = prepareLPS(pattern);

		int i = 0; // text index
		int j = 0; // pattern index

		while (i < N) {
			if (pattern[j] == text[i]) {
				j++;
				i++;
			}

			if (j == M) {
				matches.push_back(i - j);
				j = lps[j - 1];
			}

			else if (i < N && pattern[j] != text[i]) {
				if (j != 0)
					j = lps[j - 1];
				else
					i = i + 1;
			}
		}
		t.stop();

		cout << "[KMP MATCH]: Text length: " << text.size() <<
			", Pattern length: " << pattern.size() <<
			" -> Time: " << t.elapsed() << endl;
	}
	else
	{
		cout << "[KMP MATCH] Text and pattern did not load properly...";
	}
	return matches;
}


vector<int> PatternMatcher::smith_waterman()
{
	vector<int> matches;

	size_t lengthText = text.size();
	size_t lengthPattern = pattern.size();

	if (lengthText > 0 && lengthPattern > 0 && lengthPattern <= lengthText)
	{
		Timer t;
		t.start();

		size_t N = lengthText + 1;
		size_t M = lengthPattern + 1;

		double threshhold = 0.9 * lengthPattern * 2;
		cout << (int)threshhold << endl;

		int* matrix;
		int* ptr;

		matrix = new int[N * M];
		ptr = new int[N * M];

		matrix_set(matrix, M, 0, 0, 0);
		matrix_set(ptr, M, 0, 0, DONE);
		for (int i = 1; i < M; i++)
		{
			matrix_set(matrix, M, 0, i, 0);
			matrix_set(ptr, M, 0, i, LEFT);
		}
		for (int i = 1; i < N; i++)
		{
			matrix_set(matrix, M, i, 0, 0);
			matrix_set(ptr, M, i, 0, UP);
		}

		vector<int> above;
		vector<int> above_i;
		vector<int> above_j;

		for (int i = 1; i < N; i++)
		{
			for (int j = 1; j < M; j++)
			{
				int score_up = matrix_get(matrix, M, i - 1, j) + DASH_PENALTY;
				int score_left = matrix_get(matrix, M, i, j - 1) + DASH_PENALTY;
				int score_diag = matrix_get(matrix, M, i - 1, j - 1) + similarityScore(text[i - 1], pattern[j - 1]);
				double score_max = max(0, max(score_up, max(score_left, score_diag)));
				matrix_set(matrix, M, i, j, score_max);
				if (score_max == score_up) {
					matrix_set(ptr, M, i, j, UP);
				}
				else if (score_max == score_left) {
					matrix_set(ptr, M, i, j, LEFT);
				}
				else {
					matrix_set(ptr, M, i, j, DIAG);
				}

				if (score_max >= threshhold)
				{
					above.push_back(score_max);
					above_i.push_back(i);
					above_j.push_back(j);
				}
			}
		}
		for (int i = 0; i < above.size(); i++)
		{
			int len = 0;
			char* a_aligned = new char[N + M];
			char* b_aligned = new char[N + M];

			int cur_i = above_i[i];
			int cur_j = above_j[i];

			while ((cur_i != 0) && (cur_j != 0) && (matrix_get(matrix, M, cur_i, cur_j) != 0)) {
				if (matrix_get(ptr, M, cur_i, cur_j) == DIAG) {
					a_aligned[len] = text[cur_i - 1];
					b_aligned[len] = pattern[cur_j - 1];
					len++;
					cur_i--;
					cur_j--;
				}
				else if (matrix_get(ptr, M, cur_i, cur_j) == LEFT) {
					a_aligned[len] = '-';
					b_aligned[len] = pattern[cur_j - 1];
					len++;
					cur_j--;
				}
				else if (matrix_get(ptr, M, cur_i, cur_j) == UP) {
					a_aligned[len] = text[cur_i - 1];
					b_aligned[len] = '-';
					len++;
					cur_i--;
				}
			}
			matches.push_back(cur_i);
			a_aligned[len] = '\0';
			b_aligned[len] = '\0';

			for (int k = 0; k < len / 2; k++)
			{
				char tmp = a_aligned[k];
				a_aligned[k] = a_aligned[len - 1 - k];
				a_aligned[len - 1 - k] = tmp;
				tmp = b_aligned[k];
				b_aligned[k] = b_aligned[len - 1 - k];
				b_aligned[len - 1 - k] = tmp;
			}
		}

		t.stop();
		set<int> s;
		unsigned size = matches.size();
		for (unsigned i = 0; i < size; ++i) s.insert(matches[i]);
		matches.assign(s.begin(), s.end());

		cout << "[Smith-Waterman]: Text length: " << text.size() <<
			", Pattern length: " << pattern.size() <<
			" -> Time: " << t.elapsed() <<
			" Matches: " << matches.size() << endl;
	}
	else
	{
		cout << "[Smith-Waterman]: Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::naiveOpenMP()
{
	vector<int> matches;

	size_t N = text.size();
	size_t M = pattern.size();

	if (N > 0 && M > 0 && M <= N)
	{
		Timer t;
		t.start();
		int i;
#pragma omp parallel for num_threads(2) schedule(static,200) shared(matches) private(i)
		for (i = 0; i < N - M + 1; i++) {
			int j = 0;

			while (j < M && text[i + j] == pattern[j]) j++;

			if (j == M)
			{
				#pragma omp critical
				{
					matches.push_back(i);
				}
			}
		}
		t.stop();

		cout << "[NAIVE MATCH with openMP]: Text length: " << text.size() <<
			", Pattern length: " << pattern.size() <<
			" -> Time: " << t.elapsed() << endl;
	}
	else
	{
		cout << "[NAIVE MATCH with openMP]: Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::naiveParallel()
{

	size_t N = text.size();
	size_t M = pattern.size();

	vector<int> matches;

	if (N > 0 && M > 0 && M <= N)
	{
		Timer t;
		int rank;
		int size;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Status status;
		MPI_Status statuses[32];
		MPI_Request requests[32];
		int recvcnts[32];
		int localCounts = 0;
		int matchesSize = 0;
		int displs[32] = { 0 };
		vector<int> localMatches;
		if (rank == 0)
		{
			t.start();
		}

		int chunkSize = N / size;
		if (rank == 0)
		{
			int i;
			for (i = 1; i < size - 1; i++)
			{
				MPI_Send(&text[i * chunkSize - (M >> 1)], chunkSize + M - 1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
			}
			MPI_Send(&text[i * chunkSize - (M >> 1)], chunkSize + (M >> 1), MPI_CHAR, i, 0, MPI_COMM_WORLD);

			vector<char>::const_iterator first = text.begin() + 0;
			vector<char>::const_iterator last = text.begin() + chunkSize + (M >> 1);
			vector<char> data(first, last);

			localMatches = PatternMatcher::naive(data, pattern);
			recvcnts[0] = localMatches.size();

			matchesSize += recvcnts[0];
			for (int i = 1; i < size; i++)
			{
				int count;
				vector<int> buff;
				buff.resize(N - M);
				MPI_Recv(&recvcnts[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &statuses[i]);
				matchesSize += recvcnts[i];
			}
			for (int i = 1; i < size; i++)
			{
				displs[i] = recvcnts[i - 1] + displs[i - 1];
			}
		}
		else if (rank < size - 1)
		{
			vector<char> data;
			data.resize(chunkSize + M);
			MPI_Recv(&data[0], chunkSize + M - 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
			localMatches = PatternMatcher::naive(data, pattern);
			for (int i = 0; i < localMatches.size(); i++)
			{
				localMatches[i] += rank * chunkSize - (M >> 1);
			}
			localCounts = localMatches.size();
			MPI_Send(&localCounts, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		else
		{
			vector<char> data;
			data.resize(chunkSize + M / 2);
			MPI_Recv(&data[0], chunkSize + M, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
			localMatches = PatternMatcher::naive(data, pattern);
			for (int i = 0; i < localMatches.size(); i++)
			{
				localMatches[i] += rank * chunkSize - (M >> 1);
			}
			localCounts = localMatches.size();
			MPI_Send(&localCounts, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == 0)
		{
			cout << "Matches: " << matchesSize << endl;
			for (int x : localMatches)
			{
				matches.push_back(x);
			}
			for (int i = 1; i < size; i++)
			{
				if (recvcnts[i] > 0)
				{
					vector<int> recv(recvcnts[i]);
					MPI_Recv(&recv[0], recvcnts[i], MPI_INT, i, 0, MPI_COMM_WORLD, &status);
					for (int x : recv)
					{
						matches.push_back(x);
					}
				}
			}
		}
		else
		{
			if (localCounts > 0)
			{
				MPI_Send(&localMatches[0], localCounts, MPI_INT, 0, 0, MPI_COMM_WORLD);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == 0)
		{
			t.stop();
			cout << "[Naive parallel]: Text length: " << text.size() <<
				", Pattern length: " << pattern.size() <<
				" -> Time: " << t.elapsed() << endl;
		}
	}
	else
	{
		cout << "[Naive parallel] Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::naiveParallelOpenMP()
{

	size_t N = text.size();
	size_t M = pattern.size();

	vector<int> matches;

	if (N > 0 && M > 0 && M <= N)
	{
		Timer t;
		int rank;
		int size;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Status status;
		MPI_Status statuses[32];
		MPI_Request requests[32];
		int recvcnts[32];
		int localCounts = 0;
		int matchesSize = 0;
		int displs[32] = { 0 };
		vector<int> localMatches;
		if (rank == 0)
		{
			t.start();
		}

		int chunkSize = N / size;
		if (rank == 0)
		{
			int i;
			for (i = 1; i < size - 1; i++)
			{
				MPI_Send(&text[i * chunkSize - (M >> 1)], chunkSize + M - 1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
			}
			MPI_Send(&text[i * chunkSize - (M >> 1)], chunkSize + (M >> 1), MPI_CHAR, i, 0, MPI_COMM_WORLD);

			vector<char>::const_iterator first = text.begin() + 0;
			vector<char>::const_iterator last = text.begin() + chunkSize + (M >> 1);
			vector<char> data(first, last);

			localMatches = PatternMatcher::naiveOpenMP(data, pattern);
			recvcnts[0] = localMatches.size();

			matchesSize += recvcnts[0];
			for (int i = 1; i < size; i++)
			{
				int count;
				vector<int> buff;
				buff.resize(N - M);
				MPI_Recv(&recvcnts[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &statuses[i]);
				matchesSize += recvcnts[i];
			}
			for (int i = 1; i < size; i++)
			{
				displs[i] = recvcnts[i - 1] + displs[i - 1];
			}
		}
		else if (rank < size - 1)
		{
			vector<char> data;
			data.resize(chunkSize + M);
			MPI_Recv(&data[0], chunkSize + M - 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
			localMatches = PatternMatcher::naiveOpenMP(data, pattern);
			for (int i = 0; i < localMatches.size(); i++)
			{
				localMatches[i] += rank * chunkSize - (M >> 1);
			}
			localCounts = localMatches.size();
			MPI_Send(&localCounts, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		else
		{
			vector<char> data;
			data.resize(chunkSize + M / 2);
			MPI_Recv(&data[0], chunkSize + M, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
			localMatches = PatternMatcher::naiveOpenMP(data, pattern);
			for (int i = 0; i < localMatches.size(); i++)
			{
				localMatches[i] += rank * chunkSize - (M >> 1);
			}
			localCounts = localMatches.size();
			MPI_Send(&localCounts, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == 0)
		{
			cout << "Matches: " << matchesSize << endl;
			for (int x : localMatches)
			{
				matches.push_back(x);
			}
			for (int i = 1; i < size; i++)
			{
				if (recvcnts[i] > 0)
				{
					vector<int> recv(recvcnts[i]);
					MPI_Recv(&recv[0], recvcnts[i], MPI_INT, i, 0, MPI_COMM_WORLD, &status);
					for (int x : recv)
					{
						matches.push_back(x);
					}
				}
			}
		}
		else
		{
			if (localCounts > 0)
			{
				MPI_Send(&localMatches[0], localCounts, MPI_INT, 0, 0, MPI_COMM_WORLD);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == 0)
		{
			t.stop();
			cout << "[Naive parallel with OpenMP]: Text length: " << text.size() <<
				", Pattern length: " << pattern.size() <<
				" -> Time: " << t.elapsed() << endl;
		}
	}
	else
	{
		cout << "[Naive parallel with OpenMP] Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::naive(vector<char> text, vector<char> pattern)
{
	vector<int> matches;

	size_t N = text.size();
	size_t M = pattern.size();

	if (N > 0 && M > 0 && M <= N)
	{
		for (int i = 0; i < N - M + 1; i++) {
			int j = 0;

			while (j < M && text[i + j] == pattern[j]) j++;

			if (j == M)
			{
				matches.push_back(i);
			}
		}
	}
	else
	{
		cout << "[NAIVE MATCH]: Text and pattern did not load properly...";
	}
	return matches;
}

std::vector<int> PatternMatcher::naiveOpenMP(std::vector<char> text, std::vector<char> pattern)
{
	vector<int> matches;

	size_t N = text.size();
	size_t M = pattern.size();

	if (N > 0 && M > 0 && M <= N)
	{
		int i;
#pragma omp parallel for num_threads(2) schedule(static,200) shared(matches) private(i)
		for (i = 0; i < N - M + 1; i++) {
			int j = 0;

			while (j < M && text[i + j] == pattern[j]) j++;

			if (j == M)
			{
#pragma omp critical
				{
					matches.push_back(i);
				}
			}
		}
	}
	else
	{
		cout << "[NAIVE MATCH with openMP]: Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::smithWaterman(vector<char> text, vector<char> pattern)
{
	vector<int> matches;

	size_t lengthText = text.size();
	size_t lengthPattern = pattern.size();

	if (lengthText > 0 && lengthPattern > 0 && lengthPattern <= lengthText)
	{

		size_t N = lengthText + 1;
		size_t M = lengthPattern + 1;

		double threshhold = 0.9 * lengthPattern * 2;
		cout << (int)threshhold << endl;

		int* matrix;
		int* ptr;

		matrix = new int[N * M];
		ptr = new int[N * M];

		matrix_set(matrix, M, 0, 0, 0);
		matrix_set(ptr, M, 0, 0, DONE);
		for (int i = 1; i < M; i++)
		{
			matrix_set(matrix, M, 0, i, 0);
			matrix_set(ptr, M, 0, i, LEFT);
		}
		for (int i = 1; i < N; i++)
		{
			matrix_set(matrix, M, i, 0, 0);
			matrix_set(ptr, M, i, 0, UP);
		}

		vector<int> above;
		vector<int> above_i;
		vector<int> above_j;

		for (int i = 1; i < N; i++)
		{
			for (int j = 1; j < M; j++)
			{
				int score_up = matrix_get(matrix, M, i - 1, j) + DASH_PENALTY;
				int score_left = matrix_get(matrix, M, i, j - 1) + DASH_PENALTY;
				int score_diag = matrix_get(matrix, M, i - 1, j - 1) + similarityScore(text[i - 1], pattern[j - 1]);
				double score_max = max(0, max(score_up, max(score_left, score_diag)));
				matrix_set(matrix, M, i, j, score_max);
				if (score_max == score_up) {
					matrix_set(ptr, M, i, j, UP);
				}
				else if (score_max == score_left) {
					matrix_set(ptr, M, i, j, LEFT);
				}
				else {
					matrix_set(ptr, M, i, j, DIAG);
				}

				if (score_max >= threshhold)
				{
					above.push_back(score_max);
					above_i.push_back(i);
					above_j.push_back(j);
				}
			}
		}
		for (int i = 0; i < above.size(); i++)
		{
			int len = 0;
			char* a_aligned = new char[N + M];
			char* b_aligned = new char[N + M];

			int cur_i = above_i[i];
			int cur_j = above_j[i];

			while ((cur_i != 0) && (cur_j != 0) && (matrix_get(matrix, M, cur_i, cur_j) != 0)) {
				if (matrix_get(ptr, M, cur_i, cur_j) == DIAG) {
					a_aligned[len] = text[cur_i - 1];
					b_aligned[len] = pattern[cur_j - 1];
					len++;
					cur_i--;
					cur_j--;
				}
				else if (matrix_get(ptr, M, cur_i, cur_j) == LEFT) {
					a_aligned[len] = '-';
					b_aligned[len] = pattern[cur_j - 1];
					len++;
					cur_j--;
				}
				else if (matrix_get(ptr, M, cur_i, cur_j) == UP) {
					a_aligned[len] = text[cur_i - 1];
					b_aligned[len] = '-';
					len++;
					cur_i--;
				}
			}
			matches.push_back(cur_i);
			a_aligned[len] = '\0';
			b_aligned[len] = '\0';

			for (int k = 0; k < len / 2; k++)
			{
				char tmp = a_aligned[k];
				a_aligned[k] = a_aligned[len - 1 - k];
				a_aligned[len - 1 - k] = tmp;
				tmp = b_aligned[k];
				b_aligned[k] = b_aligned[len - 1 - k];
				b_aligned[len - 1 - k] = tmp;
			}
		}
		set<int> s;
		unsigned size = matches.size();
		for (unsigned i = 0; i < size; ++i) s.insert(matches[i]);
		matches.assign(s.begin(), s.end());
	}
	else
	{
		cout << "[Smith-Waterman]: Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> prepareLPS(vector<char> pattern)
{
	vector<int> lps(pattern.size());
	int len = 0;
	lps[0] = 0;

	unsigned int i = 1;
	while (i < pattern.size()) {
		if (pattern[i] == pattern[len]) {
			len++;
			lps[i] = len;
			i++;
		}
		else
		{
			if (len != 0) {
				len = lps[len - 1];
			}
			else
			{
				lps[i] = 0;
				i++;
			}
		}
	}
	return lps;
}