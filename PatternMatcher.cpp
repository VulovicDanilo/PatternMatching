#include "PatternMatcher.h"
#include "utils.h"

using namespace std;

vector<int> prepareLPS(vector<char> pattern);
vector<int> prepareBCH(vector<char> pattern);

vector<unsigned char> PatternMatcher::encode(vector<char> text)
{
	vector<unsigned char> coded;

	int N = text.size();

	int groupSize = 4;
	int numOfGroups = N / 4;
	coded.reserve(numOfGroups);

	for (int i = 0; i < N; i += 4)
	{
		char first = text[i + 0];
		char second = text[i + 1];
		char third = text[i + 2];
		char fourth = text[i + 3];

		char group = 0;

		group += getCodedValue(first);
		group <<= 2;
		group += getCodedValue(second);
		group <<= 2;
		group += getCodedValue(third);
		group <<= 2;
		group += getCodedValue(fourth);

		coded.push_back(group);
	}
	return coded;
}

vector<unsigned char> PatternMatcher::encode2(vector<char> text)
{
	vector<unsigned char> coded;
	coded.reserve(text.size() / 4);

	for (int i = 0; i < text.size(); i += 4)
	{
		char group = 0;

		group |= (text[i] >> 1) % 4;
		group <<= 2;
		group |= (text[i + 1] >> 1) % 4;
		group <<= 2;
		group |= (text[i + 2] >> 1) % 4;
		group <<= 2;
		group |= (text[i + 3] >> 1) % 4;

		coded.push_back(group);
	}
	return coded;
}

vector<unsigned char> PatternMatcher::encode3(vector<char> text)
{
	vector<unsigned char> coded;
	coded.reserve(text.size() / 4);

	for (int i = 0; i < text.size(); i += 4)
	{
		char group = 0;
		char mask = 3;
		group |= (text[i] >> 1)& mask;
		group <<= 2;
		group |= (text[i + 1] >> 1)& mask;
		group <<= 2;
		group |= (text[i + 2] >> 1)& mask;
		group <<= 2;
		group |= (text[i + 3] >> 1)& mask;

		coded.push_back(group);
	}
	return coded;
}

char PatternMatcher::getCodedValue(char c)
{
	if (c == 'A')
		return 0b00000000;
	else if (c == 'C')
		return 0b00000001;
	else if (c == 'G')
		return 0b00000010;
	else if (c == 'T')
		return 0b00000011;
	else
		return 0;
}

void PatternMatcher::loadText(std::string textFileName)
{
	Timer t;
	ifstream tfs(textFileName);
	if (!tfs.fail())
	{
		t.start();
		text.assign((istreambuf_iterator<char>(tfs)),
			(istreambuf_iterator<char>()));

		text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
		t.stop();
		cout << "Text loaded in: " << t.elapsed() << " sekundi." << endl;
	}
	else
	{
		cout << "there is a problem with one of file streams" << endl;
	}
}

void PatternMatcher::loadPattern(std::string patternFileName)
{
	Timer t;
	ifstream pfs(patternFileName);
	if (!pfs.fail())
	{
		t.start();
		pattern.assign((istreambuf_iterator<char>(pfs)),
			(istreambuf_iterator<char>()));

		pattern.erase(std::remove(pattern.begin(), pattern.end(), '\n'), pattern.end());
		pattern.erase(std::remove(pattern.begin(), pattern.end(), '\r'), pattern.end());
		t.stop();
		cout << "Pattern loaded in: " << t.elapsed() << " sekundi." << endl;
	}
	else
	{
		cout << "there is a problem with one of file streams" << endl;
	}
}

void PatternMatcher::loadTextChunk(string textFileName)
{
	Timer t;
	t.start();
	ifstream ifs(textFileName);
	char buff[FASTA_FILE_ROW_SIZE];
	while (ifs.read(buff, FASTA_FILE_ROW_SIZE))
	{
		text.insert(text.end(), buff, buff + 64);
	}
	t.stop();
	cout << "Text loaded in " << t.elapsed() << " seconds." << endl;
}

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

		cout << "[NAIVE]: Text length: " << text.size() <<
			", Pattern length: " << pattern.size() <<
			" -> Time: " << t.elapsed() << endl;
	}
	else
	{
		cout << "[NAIVE]: Text and pattern did not load properly...";
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
		Timer t;
		t.start();

		vector<int> badChars = prepareBCH(pattern);

		int shift = 0;
		while (shift <= (N - M))
		{
			int j = M - 1;
			while (j >= 0 && pattern[j] == text[shift + j])
				j--;
			if (j < 0)
			{
				matches.push_back(shift);
				shift += (shift + M < N) ? M - badChars[text[shift + M]] : 1;
			}
			else
				shift += max(1, j - badChars[text[shift + j]]);
		}
		t.stop();

		cout << "[BOYER-MOORE]: Text length: " << text.size() <<
			", Pattern length: " << pattern.size() <<
			" -> Time: " << t.elapsed() << endl;
	}
	else
	{
		cout << "[BOYER-MOORE] Text and pattern did not load properly...";
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

		int i = 0;
		int j = 0;

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

		cout << "[KMP]: Text length: " << text.size() <<
			", Pattern length: " << pattern.size() <<
			" -> Time: " << t.elapsed() << endl;
	}
	else
	{
		cout << "[KMP] Text and pattern did not load properly...";
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
				int score_up = matrix_get(matrix, M, i - 1, j) + GAP_PENALTY;
				int score_left = matrix_get(matrix, M, i, j - 1) + GAP_PENALTY;
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

vector<int> PatternMatcher::needleman_wunsch()
{
	// TODO
	vector<int> matches;


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

		cout << "[NAIVE with openMP]: Text length: " << text.size() <<
			", Pattern length: " << pattern.size() <<
			" -> Time: " << t.elapsed() << endl;
	}
	else
	{
		cout << "[NAIVE with openMP]: Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::naiveParallel()
{
	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	size_t N = text.size();
	size_t M = pattern.size();

	int chunkSize = N / size;

	MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&M, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&chunkSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank)
	{
		pattern.resize(M);
	}
	MPI_Bcast(&pattern[0], pattern.size(), MPI_CHAR, 0, MPI_COMM_WORLD);
	vector<int> matches;

	if (N > 0 && M > 0 && M <= N)
	{
		Timer t;
		MPI_Status status;
		int* counts = new int[size];
		int* displs = new int[size];
		int* alldata = new int[4];
		int localCounts = 0;
		vector<int> localMatches;
		for (int i = 0; i < size; i++)
		{
			counts[i] = displs[i] = 0;
		}
		if (rank == 0)
		{
			t.start();
		}
		if (rank == 0)
		{
			int i;
			for (i = 1; i < size - 1; i++)
			{
				MPI_Send(&text[i * chunkSize - (M >> 1)], chunkSize + M, MPI_CHAR, i, 0, MPI_COMM_WORLD);
			}
			MPI_Send(&text[i * chunkSize - (M >> 1)], chunkSize + (M >> 1), MPI_CHAR, i, 0, MPI_COMM_WORLD);

			vector<char> data(text.begin(), text.begin() + chunkSize + (M >> 1));
			localMatches = PatternMatcher::naive(data, pattern);
		}
		else
		{
			vector<char> buff;
			buff.resize(chunkSize + M);
			MPI_Recv(buff.data(), chunkSize + M, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
			localMatches = PatternMatcher::naive(buff, pattern);
		}
		for (int i = 0; i < localMatches.size(); i++)
		{
			if (rank)
				localMatches[i] += rank * chunkSize - (M >> 1);
		}
		localCounts = localMatches.size();
		MPI_Gather(&localCounts, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			for (int i = 0; i < size; i++)
			{
				displs[i] = (i > 0) ? (displs[i - 1] + counts[i - 1]) : 0;
			}
			alldata = new int[counts[size - 1] + displs[size - 1]];
			for (int k = 0; k < counts[size - 1] + displs[size - 1]; k++)
				alldata[k] = 0;
		}
		MPI_Gatherv(localMatches.data(), localCounts, MPI_INT, alldata, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			matches.assign(alldata, alldata + counts[size - 1] + displs[size - 1]);
			t.stop();
			cout << "[Naive parallel]: Text length: " << text.size() <<
				", Pattern length: " << pattern.size() <<
				" -> Time: " << t.elapsed() << endl;
		}
		delete[] counts;
		delete[] displs;
		delete[] alldata;
	}
	else
	{
		cout << "[Naive parallel] Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::naiveParallelOpenMP()
{
	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	size_t N = text.size();
	size_t M = pattern.size();

	int chunkSize = N / size;

	MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&M, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&chunkSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank)
	{
		pattern.resize(M);
	}
	MPI_Bcast(&pattern[0], pattern.size(), MPI_CHAR, 0, MPI_COMM_WORLD);
	vector<int> matches;

	if (N > 0 && M > 0 && M <= N)
	{
		Timer t;
		MPI_Status status;
		int* counts = new int[size];
		int* displs = new int[size];
		int* alldata = new int[4];
		int localCounts = 0;
		vector<int> localMatches;
		if (rank == 0)
		{
			t.start();
		}
		if (rank == 0)
		{
			int i;
			for (i = 1; i < size - 1; i++)
			{
				MPI_Send(&text[i * chunkSize - (M >> 1)], chunkSize + M - 1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
			}
			MPI_Send(&text[i * chunkSize - (M >> 1)], chunkSize + (M >> 1), MPI_CHAR, i, 0, MPI_COMM_WORLD);

			vector<char> data(text.begin(), text.begin() + chunkSize + (M >> 1));
			localMatches = PatternMatcher::naiveOpenMP(data, pattern);
		}
		else
		{
			vector<char> buff;
			buff.resize(chunkSize + M);
			MPI_Recv(buff.data(), chunkSize + M, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
			localMatches = PatternMatcher::naiveOpenMP(buff, pattern);
		}
		for (int i = 0; i < localMatches.size(); i++)
		{
			if (rank)
				localMatches[i] += rank * chunkSize - (M >> 1);
		}
		localCounts = localMatches.size();
		MPI_Gather(&localCounts, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			for (int i = 0; i < size; i++)
			{
				displs[i] = (i > 0) ? (displs[i - 1] + counts[i - 1]) : 0;
			}
			alldata = new int[counts[size - 1] + displs[size - 1]];
			for (int k = 0; k < counts[size - 1] + displs[size - 1]; k++)
				alldata[k] = 0;
		}
		MPI_Gatherv(localMatches.data(), localCounts, MPI_INT, alldata, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			matches.assign(alldata, alldata + counts[size - 1] + displs[size - 1]);
			t.stop();
			cout << "[Naive parallel with OpenMP]: Text length: " << text.size() <<
				", Pattern length: " << pattern.size() <<
				" -> Time: " << t.elapsed() << endl;
		}
		delete[] counts;
		delete[] displs;
		delete[] alldata;
	}
	else
	{
		cout << "[Naive parallel with OpenMP] Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::boyer_mooreParallel()
{
	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	size_t N = text.size();
	size_t M = pattern.size();

	int chunkSize = N / size;

	MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&M, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&chunkSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank)
	{
		pattern.resize(M);
	}
	MPI_Bcast(&pattern[0], pattern.size(), MPI_CHAR, 0, MPI_COMM_WORLD);
	vector<int> matches;

	if (N > 0 && M > 0 && M <= N)
	{
		Timer t;
		MPI_Status status;
		int* counts = new int[size];
		int* displs = new int[size];
		int* alldata = new int[4];
		for (int i = 0; i < size; i++)
		{
			counts[i] = displs[i] = 0;
		}
		int localCounts = 0;
		vector<int> localMatches;
		if (rank == 0)
		{
			t.start();
		}
		if (rank == 0)
		{
			int i;
			for (i = 1; i < size - 1; i++)
			{
				MPI_Send(&text[i * chunkSize - (M >> 1)], chunkSize + M, MPI_CHAR, i, 0, MPI_COMM_WORLD);
			}
			MPI_Send(&text[i * chunkSize - (M >> 1)], chunkSize + (M >> 1), MPI_CHAR, i, 0, MPI_COMM_WORLD);

			vector<char> data(text.begin(), text.begin() + chunkSize + (M >> 1));
			localMatches = PatternMatcher::boyer_moore(data, pattern);
		}
		else
		{
			vector<char> buff;
			buff.resize(chunkSize + M);
			MPI_Recv(buff.data(), chunkSize + M, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
			localMatches = PatternMatcher::boyer_moore(buff, pattern);
		}
		for (int i = 0; i < localMatches.size(); i++)
		{
			if (rank)
				localMatches[i] += rank * chunkSize - (M >> 1);
		}
		localCounts = localMatches.size();
		MPI_Gather(&localCounts, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			for (int i = 0; i < size; i++)
			{
				displs[i] = (i > 0) ? (displs[i - 1] + counts[i - 1]) : 0;
			}
			alldata = new int[counts[size - 1] + displs[size - 1]];
			for (int k = 0; k < counts[size - 1] + displs[size - 1]; k++)
				alldata[k] = 0;
		}
		MPI_Gatherv(localMatches.data(), localCounts, MPI_INT, alldata, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			matches.assign(alldata, alldata + counts[size - 1] + displs[size - 1]);
			t.stop();
			cout << "[BOYER-MOORE parallel]: Text length: " << text.size() <<
				", Pattern length: " << pattern.size() <<
				" -> Time: " << t.elapsed() << endl;
		}
		delete[] counts;
		delete[] displs;
		delete[] alldata;
	}
	else
	{
		cout << "[BOYER-MOORE parallel] Text and pattern did not load properly...";
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
		cout << "[NAIVE]: Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::boyer_moore(vector<char> text, vector<char> pattern)
{
	vector<int> matches;

	size_t N = text.size();
	size_t M = pattern.size();

	if (N > 0 && M > 0 && M <= N)
	{

		vector<int> badChars = prepareBCH(pattern);

		int shift = 0;
		while (shift <= (N - M))
		{
			int j = M - 1;
			while (j >= 0 && pattern[j] == text[shift + j])
				j--;
			if (j < 0)
			{
				matches.push_back(shift);
				shift += (shift + M < N) ? M - badChars[text[shift + M]] : 1;
			}

			else
				shift += max(1, j - badChars[text[shift + j]]);
		}
	}
	else
	{
		cout << "[BOYER-MOORE] Text and pattern did not load properly...";
	}

	return matches;
}

vector<int> PatternMatcher::naiveOpenMP(vector<char> text, vector<char> pattern)
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
		cout << "[NAIVE with openMP]: Text and pattern did not load properly...";
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
				int score_up = matrix_get(matrix, M, i - 1, j) + GAP_PENALTY;
				int score_left = matrix_get(matrix, M, i, j - 1) + GAP_PENALTY;
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

vector<int> PatternMatcher::coded_naive()
{
	vector<int> matches;

	vector<unsigned char> codedText = encode3(text);
	vector<unsigned char> codedPattern = encode3(pattern);

	size_t N = codedText.size();
	size_t M = codedPattern.size();

	if (N > 0 && M > 0 && M <= N)
	{
		CodedPattern patterns[4];
		Timer t;
		t.start();
		CodedPattern pat1(codedPattern, 0);
		CodedPattern pat2(codedPattern, 1);
		CodedPattern pat3(codedPattern, 2);
		CodedPattern pat4(codedPattern, 3);

		patterns[0] = pat1;
		patterns[1] = pat2;
		patterns[2] = pat3;
		patterns[3] = pat4;
		for (int i = 0; i < N - M; i++)
		{
			int j;
			for (auto pattern : patterns)
			{
				j = 1;
				while (j < pattern.size() - 1 && codedText[i + j] == pattern[j]) j++;
				if (j == pattern.size() - 1)
				{
					if (pattern.compareFirstByte(codedText[i]) && pattern.compareLastByte(codedText[i + j]))
						matches.push_back(i * 4 + pattern.prefix);
				}
			}
		}

		t.stop();
		cout << "[NAIVE - CODED]: Text length: " << text.size() <<
			", Pattern length: " << pattern.size() <<
			" -> Time: " << t.elapsed() << endl;
	}
	else
	{
		cout << "[NAIVE - CODED]: Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::coded_boyer_moore()
{
	vector<int> matches;

	vector<unsigned char> codedText = encode3(text);
	vector<unsigned char> codedPattern = encode3(pattern);

	size_t N = codedText.size();
	size_t M = codedPattern.size();

	if (N > 0 && M > 0 && M <= N)
	{
		CodedPattern patterns[4];
		Timer t;
		t.start();
		CodedPattern pat1(codedPattern, 0);
		CodedPattern pat2(codedPattern, 1);
		CodedPattern pat3(codedPattern, 2);
		CodedPattern pat4(codedPattern, 3);

		patterns[0] = pat1;
		patterns[1] = pat2;
		patterns[2] = pat3;
		patterns[3] = pat4;

		for (int k = 0; k < 4; k++)
		{
			int shift = 0;
			int localM = patterns[k].size();
			while (shift <= (N - localM))
			{
				int j = localM - 2;
				while (j > 0 && patterns[k].content[j] == codedText[shift + j])
					j--;
				if (j == 0)
				{
					if (patterns[k].compareFirstByte(codedText[shift]) && patterns[k].compareLastByte(codedText[shift + localM - 1]))
						matches.push_back((shift << 2) + patterns[k].prefix);
					shift += (shift + localM < N) ? localM - patterns[k].badChars[codedText[shift + localM]] : 1;
				}
				else
				{
					shift += max(1, j - 1 - patterns[k].badChars[codedText[shift + j]]);
				}
			}
		}

		t.stop();
		cout << "[BOYER-MOORE - CODED]: Text length: " << text.size() <<
			", Pattern length: " << pattern.size() <<
			" -> Time: " << t.elapsed() << endl;
	}
	else
	{
		cout << "[BOYER-MOORE - CODED]: Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::coded_naiveOpenMP()
{
	vector<int> matches;

	vector<unsigned char> codedText = encode3(text);
	vector<unsigned char> codedPattern = encode3(pattern);

	size_t N = codedText.size();
	size_t M = codedPattern.size();

	if (N > 0 && M > 0 && M <= N)
	{
		CodedPattern patterns[4];
		int i;

		Timer t;
		t.start();
		CodedPattern pat1(codedPattern, 0);
		CodedPattern pat2(codedPattern, 1);
		CodedPattern pat3(codedPattern, 2);
		CodedPattern pat4(codedPattern, 3);

		patterns[0] = pat1;
		patterns[1] = pat2;
		patterns[2] = pat3;
		patterns[3] = pat4;
#pragma omp parallel for num_threads(2) schedule(static,100) shared(matches) private(i)
		for (i = 0; i < N - M; i++)
		{
			int j;
			for (auto pattern : patterns)
			{
				j = 1;
				while (j < pattern.size() - 1 && codedText[i + j] == pattern[j]) j++;
				if (j == pattern.size() - 1)
				{
					if (pattern.compareFirstByte(codedText[i]) && pattern.compareLastByte(codedText[i + j]))
#pragma omp critical
					{
						matches.push_back(i * 4 + pattern.prefix);
					}
				}
			}
		}
		t.stop();

		cout << "[NAIVE with openMP - CODED]: Text length: " << text.size() <<
			", Pattern length: " << pattern.size() <<
			" -> Time: " << t.elapsed() << endl;
	}
	else
	{
		cout << "[NAIVE with openMP - CODED]: Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::coded_naiveParallel()
{
	vector<unsigned char> codedText = encode3(text);
	vector<unsigned char> codedPattern = encode3(pattern);

	size_t N = codedText.size();
	size_t M = codedPattern.size();

	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int chunkSize = N / size;

	MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&M, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&chunkSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank)
	{
		codedPattern.resize(M);
	}
	MPI_Bcast(&codedPattern[0], M, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
	vector<int> matches;

	if (N > 0 && M > 0 && M <= N)
	{
		Timer t;
		MPI_Status status;
		int* counts = new int[size];
		int* displs = new int[size];
		int* alldata = new int[4];
		for (int i = 0; i < size; i++)
		{
			counts[i] = displs[i] = 0;
		}
		int localCounts = 0;
		vector<int> localMatches;
		if (rank == 0)
		{
			t.start();
		}
		if (rank == 0)
		{
			int i;
			for (i = 1; i < size - 1; i++)
			{
				MPI_Send(&codedText[i * chunkSize - (M >> 1)], chunkSize + M, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);
			}
			MPI_Send(&codedText[i * chunkSize - (M >> 1)], chunkSize + (M >> 1), MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);

			vector<unsigned char> data(codedText.begin(), codedText.begin() + chunkSize + (M >> 1));
			localMatches = PatternMatcher::coded_naive(data, codedPattern);
		}
		else
		{
			vector<unsigned char> buff;
			buff.resize(chunkSize + M);
			MPI_Recv(buff.data(), chunkSize + M, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &status);
			localMatches = PatternMatcher::coded_naive(buff, codedPattern);
		}
		for (int i = 0; i < localMatches.size(); i++)
		{
			if (rank)
				localMatches[i] += 4 * (rank * chunkSize - (M >> 1));
		}
		localCounts = localMatches.size();
		MPI_Gather(&localCounts, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			for (int i = 0; i < size; i++)
			{
				displs[i] = (i > 0) ? (displs[i - 1] + counts[i - 1]) : 0;
			}
			alldata = new int[counts[size - 1] + displs[size - 1]];
			for (int k = 0; k < counts[size - 1] + displs[size - 1]; k++)
				alldata[k] = 0;
		}
		MPI_Gatherv(localMatches.data(), localCounts, MPI_INT, alldata, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			matches.assign(alldata, alldata + counts[size - 1] + displs[size - 1]);
			t.stop();
			cout << "[Naive parallel - CODED]: Text length: " << text.size() <<
				", Pattern length: " << pattern.size() <<
				" -> Time: " << t.elapsed() << endl;
		}
		delete[] counts;
		delete[] displs;
		delete[] alldata;
	}
	else
	{
		cout << "[Naive parallel - CODED] Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::coded_naiveParallelOpenMP()
{
	vector<unsigned char> codedText = encode3(text);
	vector<unsigned char> codedPattern = encode3(pattern);

	size_t N = codedText.size();
	size_t M = codedPattern.size();

	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int chunkSize = N / size;

	MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&M, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&chunkSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank)
	{
		codedPattern.resize(M);
	}
	MPI_Bcast(&codedPattern[0], M, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
	vector<int> matches;

	if (N > 0 && M > 0 && M <= N)
	{
		Timer t;
		MPI_Status status;
		int* counts = new int[size];
		int* displs = new int[size];
		int* alldata = new int[4];
		for (int i = 0; i < size; i++)
		{
			counts[i] = displs[i] = 0;
		}
		int localCounts = 0;
		vector<int> localMatches;
		if (rank == 0)
		{
			t.start();
		}

		if (rank == 0)
		{
			int i;
			for (i = 1; i < size - 1; i++)
			{
				MPI_Send(&codedText[i * chunkSize - (M >> 1)], chunkSize + M, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);
			}
			MPI_Send(&codedText[i * chunkSize - (M >> 1)], chunkSize + (M >> 1), MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);

			vector<unsigned char> data(codedText.begin(), codedText.begin() + chunkSize + (M >> 1));
			localMatches = PatternMatcher::coded_naiveOpenMP(data, codedPattern);
		}
		else
		{
			vector<unsigned char> buff;
			buff.resize(chunkSize + M);
			MPI_Recv(buff.data(), chunkSize + M, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &status);
			localMatches = PatternMatcher::coded_naiveOpenMP(buff, codedPattern);
		}
		for (int i = 0; i < localMatches.size(); i++)
		{
			if (rank)
				localMatches[i] += 4 * (rank * chunkSize - (M >> 1));
		}
		localCounts = localMatches.size();
		MPI_Gather(&localCounts, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			for (int i = 0; i < size; i++)
			{
				displs[i] = (i > 0) ? (displs[i - 1] + counts[i - 1]) : 0;
			}
			alldata = new int[counts[size - 1] + displs[size - 1]];
			for (int k = 0; k < counts[size - 1] + displs[size - 1]; k++)
				alldata[k] = 0;
		}
		MPI_Gatherv(localMatches.data(), localCounts, MPI_INT, alldata, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			matches.assign(alldata, alldata + counts[size - 1] + displs[size - 1]);
			t.stop();
			cout << "[Naive parallel with OpenMP - CODED]: Text length: " << text.size() <<
				", Pattern length: " << pattern.size() <<
				" -> Time: " << t.elapsed() << endl;
		}
		delete[] counts;
		delete[] displs;
		delete[] alldata;
	}
	else
	{
		cout << "[Naive parallel with OpenMP - CODED] Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::coded_boyer_mooreOpenMP()
{
	vector<int> matches;

	vector<unsigned char> codedText = encode3(text);
	vector<unsigned char> codedPattern = encode3(pattern);

	size_t N = codedText.size();
	size_t M = codedPattern.size();

	if (N > 0 && M > 0 && M <= N)
	{
		CodedPattern patterns[4];
		Timer t;
		t.start();
		CodedPattern pat1(codedPattern, 0);
		CodedPattern pat2(codedPattern, 1);
		CodedPattern pat3(codedPattern, 2);
		CodedPattern pat4(codedPattern, 3);

		patterns[0] = pat1;
		patterns[1] = pat2;
		patterns[2] = pat3;
		patterns[3] = pat4;
		int k = 0;
#pragma omp parallel for num_threads(4) schedule(static,1) shared(matches) private(k)
		for (k = 0; k < 4; k++)
		{
			int shift = 0;
			int localM = patterns[k].size();
			while (shift <= (N - localM))
			{
				int j = localM - 2;
				while (j > 0 && patterns[k].content[j] == codedText[shift + j])
					j--;
				if (j == 0)
				{
					if (patterns[k].compareFirstByte(codedText[shift]) && patterns[k].compareLastByte(codedText[shift + localM - 1]))
#pragma omp critical
					{
						matches.push_back((shift << 2) + patterns[k].prefix);
					}
					shift += (shift + localM < N) ? localM - patterns[k].badChars[codedText[shift + localM]] : 1;
				}
				else
				{
					shift += max(1, j - 1 - patterns[k].badChars[codedText[shift + j]]);
				}
			}
		}

		t.stop();
		cout << "[BOYER-MOORE with OpenMP - CODED]: Text length: " << text.size() <<
			", Pattern length: " << pattern.size() <<
			" -> Time: " << t.elapsed() << endl;
	}
	else
	{
		cout << "[BOYER-MOORE with OpenMP - CODED]: Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::coded_boyer_mooreParallel()
{
	Timer e;
	e.start();
	vector<unsigned char> codedText = encode3(text);
	vector<unsigned char> codedPattern = encode3(pattern);
	e.stop();

	size_t N = codedText.size();
	size_t M = codedPattern.size();

	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int chunkSize = N / size;

	MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&M, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&chunkSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank)
	{
		codedPattern.resize(M);
	}
	MPI_Bcast(&codedPattern[0], M, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
	vector<int> matches;

	if (N > 0 && M > 0 && M <= N)
	{
		Timer t;
		MPI_Status status;
		int* counts = new int[size];
		int* displs = new int[size];
		int* alldata = new int[4];
		for (int i = 0; i < size; i++)
		{
			counts[i] = displs[i] = 0;
		}
		int localCounts = 0;
		vector<int> localMatches;
		if (rank == 0)
		{
			t.start();
		}
		if (rank == 0)
		{
			int i;
			for (i = 1; i < size - 1; i++)
			{
				MPI_Send(&codedText[i * chunkSize - (M >> 1)], chunkSize + M, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);
			}
			MPI_Send(&codedText[i * chunkSize - (M >> 1)], chunkSize + (M >> 1), MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);

			vector<unsigned char> data(codedText.begin(), codedText.begin() + chunkSize + (M >> 1));
			localMatches = PatternMatcher::coded_boyer_moore(data, codedPattern);
		}
		else
		{
			vector<unsigned char> buff;
			buff.resize(chunkSize + M);
			MPI_Recv(buff.data(), chunkSize + M, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &status);
			localMatches = PatternMatcher::coded_boyer_moore(buff, codedPattern);
		}
		for (int i = 0; i < localMatches.size(); i++)
		{
			if (rank)
				localMatches[i] += 4 * (rank * chunkSize - (M >> 1));
		}
		localCounts = localMatches.size();
		MPI_Gather(&localCounts, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			for (int i = 0; i < size; i++)
			{
				displs[i] = (i > 0) ? (displs[i - 1] + counts[i - 1]) : 0;
			}
			alldata = new int[counts[size - 1] + displs[size - 1]];
			for (int k = 0; k < counts[size - 1] + displs[size - 1]; k++)
				alldata[k] = 0;
		}
		MPI_Gatherv(localMatches.data(), localCounts, MPI_INT, alldata, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			matches.assign(alldata, alldata + counts[size - 1] + displs[size - 1]);
			t.stop();
			cout << "[BOYER-MOORE parallel - CODED]: Text length: " << text.size() <<
				", Pattern length: " << pattern.size() <<
				" -> Time: " << t.elapsed() << " + Encoding time: " << e.elapsed() << endl;
		}
		delete[] counts;
		delete[] displs;
		delete[] alldata;
	}
	else
	{
		cout << "[BOYER-MOORE parallel - CODED] Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::coded_boyer_mooreParallelOpenMP()
{
	Timer e;
	e.start();
	vector<unsigned char> codedText = encode3(text);
	vector<unsigned char> codedPattern = encode3(pattern);
	e.stop();

	size_t N = codedText.size();
	size_t M = codedPattern.size();

	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int chunkSize = N / size;

	MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&M, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&chunkSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank)
	{
		codedPattern.resize(M);
	}
	MPI_Bcast(&codedPattern[0], M, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
	vector<int> matches;

	if (N > 0 && M > 0 && M <= N)
	{
		Timer t;
		MPI_Status status;
		int* counts = new int[size];
		int* displs = new int[size];
		int* alldata = new int[4];
		for (int i = 0; i < size; i++)
		{
			counts[i] = displs[i] = 0;
		}
		int localCounts = 0;
		vector<int> localMatches;
		if (rank == 0)
		{
			t.start();
		}
		if (rank == 0)
		{
			int i;
			for (i = 1; i < size - 1; i++)
			{
				MPI_Send(&codedText[i * chunkSize - (M >> 1)], chunkSize + M, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);
			}
			MPI_Send(&codedText[i * chunkSize - (M >> 1)], chunkSize + (M >> 1), MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);

			vector<unsigned char> data(codedText.begin(), codedText.begin() + chunkSize + (M >> 1));
			localMatches = PatternMatcher::coded_boyer_mooreOpenMP(data, codedPattern);
		}
		else
		{
			vector<unsigned char> buff;
			buff.resize(chunkSize + M);
			MPI_Recv(buff.data(), chunkSize + M, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &status);
			localMatches = PatternMatcher::coded_boyer_mooreOpenMP(buff, codedPattern);
		}
		for (int i = 0; i < localMatches.size(); i++)
		{
			if (rank)
				localMatches[i] += 4 * (rank * chunkSize - (M >> 1));
		}
		localCounts = localMatches.size();
		MPI_Gather(&localCounts, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			for (int i = 0; i < size; i++)
			{
				displs[i] = (i > 0) ? (displs[i - 1] + counts[i - 1]) : 0;
			}
			alldata = new int[counts[size - 1] + displs[size - 1]];
			for (int k = 0; k < counts[size - 1] + displs[size - 1]; k++)
				alldata[k] = 0;
		}
		MPI_Gatherv(localMatches.data(), localCounts, MPI_INT, alldata, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			matches.assign(alldata, alldata + counts[size - 1] + displs[size - 1]);
			t.stop();
			cout << "[BOYER-MORE parallel with OpenMP - CODED]: Text length: " << text.size() <<
				", Pattern length: " << pattern.size() <<
				" -> Time: " << t.elapsed() <<  " + Encoding time: " << e.elapsed()  << endl;
		}
		delete[] counts;
		delete[] displs;
		delete[] alldata;
	}
	else
	{
		cout << "[BOYER-MORE parallel with OpenMP - CODED]] Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::coded_naive(vector<unsigned char> codedText, vector<unsigned char> codedPattern)
{
	vector<int> matches;

	size_t N = codedText.size();
	size_t M = codedPattern.size();

	if (N > 0 && M > 0 && M <= N)
	{
		CodedPattern patterns[4];
		CodedPattern pat1(codedPattern, 0);
		CodedPattern pat2(codedPattern, 1);
		CodedPattern pat3(codedPattern, 2);
		CodedPattern pat4(codedPattern, 3);

		patterns[0] = pat1;
		patterns[1] = pat2;
		patterns[2] = pat3;
		patterns[3] = pat4;
		for (int i = 0; i < N - M; i++)
		{
			int j;
			for (auto pattern : patterns)
			{
				j = 1;
				while (j < pattern.size() - 1 && codedText[i + j] == pattern[j]) j++;
				if (j == pattern.size() - 1)
				{
					if (pattern.compareFirstByte(codedText[i]) && pattern.compareLastByte(codedText[i + j]))
					{
						matches.push_back(4 * i + pattern.prefix);
					}
				}
			}
		}
	}
	else
	{
		cout << "[NAIVE - CODED]: Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::coded_naiveOpenMP(vector<unsigned char> codedText, vector<unsigned char> codedPattern)
{
	vector<int> matches;

	size_t N = codedText.size();
	size_t M = codedPattern.size();

	if (N > 0 && M > 0 && M <= N)
	{
		CodedPattern patterns[4];
		int i;

		CodedPattern pat1(codedPattern, 0);
		CodedPattern pat2(codedPattern, 1);
		CodedPattern pat3(codedPattern, 2);
		CodedPattern pat4(codedPattern, 3);

		patterns[0] = pat1;
		patterns[1] = pat2;
		patterns[2] = pat3;
		patterns[3] = pat4;
#pragma omp parallel for num_threads(2) schedule(static,100) shared(matches) private(i)
		for (i = 0; i < N - M; i++)
		{
			int j;
			for (auto pattern : patterns)
			{
				j = 1;
				while (j < pattern.size() - 1 && codedText[i + j] == pattern[j]) j++;
				if (j == pattern.size() - 1)
				{
					if (pattern.compareFirstByte(codedText[i]) && pattern.compareLastByte(codedText[i + j]))
#pragma omp critical
					{
						matches.push_back(i * 4 + pattern.prefix);
					}
				}
			}
		}
	}
	else
	{
		cout << "[NAIVE with openMP - CODED]: Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::coded_boyer_moore(vector<unsigned char> codedText, vector<unsigned char> codedPattern)
{
	vector<int> matches;

	size_t N = codedText.size();
	size_t M = codedPattern.size();

	if (N > 0 && M > 0 && M <= N)
	{
		CodedPattern patterns[4];
		Timer t;
		t.start();
		CodedPattern pat1(codedPattern, 0);
		CodedPattern pat2(codedPattern, 1);
		CodedPattern pat3(codedPattern, 2);
		CodedPattern pat4(codedPattern, 3);

		patterns[0] = pat1;
		patterns[1] = pat2;
		patterns[2] = pat3;
		patterns[3] = pat4;

		for (int k = 0; k < 4; k++)
		{
			int shift = 0;
			int localM = patterns[k].size();
			while (shift <= (N - localM))
			{
				int j = localM - 2;
				while (j > 0 && patterns[k].content[j] == codedText[shift + j])
					j--;
				if (j == 0)
				{
					if (patterns[k].compareFirstByte(codedText[shift]) && patterns[k].compareLastByte(codedText[shift + localM - 1]))
						matches.push_back((shift << 2) + patterns[k].prefix);
					shift += (shift + localM < N) ? localM - patterns[k].badChars[codedText[shift + localM]] : 1;
				}
				else
				{
					shift += max(1, j - 1 - patterns[k].badChars[codedText[shift + j]]);
				}
			}
		}
	}
	else
	{
		cout << "[BOYER-MOORE - CODED]: Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> PatternMatcher::coded_boyer_mooreOpenMP(vector<unsigned char> codedText, vector<unsigned char> codedPattern)
{
	vector<int> matches;

	size_t N = codedText.size();
	size_t M = codedPattern.size();

	if (N > 0 && M > 0 && M <= N)
	{
		CodedPattern patterns[4];
		CodedPattern pat1(codedPattern, 0);
		CodedPattern pat2(codedPattern, 1);
		CodedPattern pat3(codedPattern, 2);
		CodedPattern pat4(codedPattern, 3);

		patterns[0] = pat1;
		patterns[1] = pat2;
		patterns[2] = pat3;
		patterns[3] = pat4;
		int k = 0;
#pragma omp parallel for num_threads(4) schedule(static,1) shared(matches) private(k)
		for (k = 0; k < 4; k++)
		{
			int shift = 0;
			int localM = patterns[k].size();
			while (shift <= (N - localM))
			{
				int j = localM - 2;
				while (j > 0 && patterns[k].content[j] == codedText[shift + j])
					j--;
				if (j == 0)
				{
					if (patterns[k].compareFirstByte(codedText[shift]) && patterns[k].compareLastByte(codedText[shift + localM - 1]))
#pragma omp critical
					{
						matches.push_back((shift << 2) + patterns[k].prefix);
					}
					shift += (shift + localM < N) ? localM - patterns[k].badChars[codedText[shift + localM]] : 1;
				}
				else
				{
					shift += max(1, j - 1 - patterns[k].badChars[codedText[shift + j]]);
				}
			}
		}
	}
	else
	{
		cout << "[BOYER-MOORE - CODED]: Text and pattern did not load properly...";
	}
	return matches;
}

vector<int> prepareBCH(vector<char> pattern)
{
	int i;
	vector<int> badChars;
	for (i = 0; i < NO_OF_CHARS; i++)
		badChars.push_back(-1);
	for (i = 0; i < pattern.size(); i++)
		badChars[(int)pattern[i]] = i;
	return badChars;
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