
#include <vector>
#include <bitset>

class CodedPattern
{
private:

	unsigned char getMask()
	{
		if (prefix == 0)
			return 0xFF;
		else if (prefix == 1)
		{
			return 0x03;
		}
		else if (prefix == 2)
		{
			return 0x0F;
		}
		else if (prefix == 3)
		{
			return 0x3F;
		}
	}

	unsigned char getFirstByteMask()
	{
		if (prefix == 0)
			return 0xFF;
		else if (prefix == 1)
		{
			return 0x3F;
		}
		else if (prefix == 2)
		{
			return 0x0F;
		}
		else if (prefix == 3)
		{
			return 0x03;
		}
	}

public:
	int prefix;
	int suffix;
	std::vector<unsigned char> content;
	unsigned char firstByte;
	unsigned char lastByte;

	inline size_t size() { return content.size(); }

	unsigned char& operator[] (int index)
	{
		return content[index];
	}
	CodedPattern()
	{
		prefix = suffix = 0;
	}
	CodedPattern(std::vector<unsigned char> pattern, int Prefix)
	{
		content.reserve(pattern.size() + 1);
		prefix = Prefix;
		suffix = 4 - Prefix;

		for (char c : pattern)
		{
			content.push_back(c);
		}
		content.push_back(0);

		if (prefix != 0)
		{
			unsigned char mask = 0;
			unsigned char oppositeMask = 0;
			mask = getMask();
			unsigned char bits1 = 0;
			unsigned char bits2 = 0;
			int i;
			for (i = 0; i < content.size(); i++)
			{
				bits1 = content[i] & mask;
				content[i] >>= (prefix * 2);
				bits2 <<= (suffix * 2);
				content[i] |= bits2;
				bits2 = bits1;
			}
		}

		firstByte = content.front();
		lastByte = content.back();
	}

	bool compareFirstByte(unsigned char b)
	{
		unsigned char mask = getFirstByteMask();
		unsigned char bits = b & mask;
		return firstByte == bits;
	}
	bool compareLastByte(unsigned char b)
	{
		unsigned char mask = getFirstByteMask();
		mask = ~mask;
		unsigned char bits = b & mask;
		return lastByte == bits;
	}
};