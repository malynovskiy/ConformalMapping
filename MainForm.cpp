#include "MainForm.h"

template <typename T>
inline T square(T a)
{
	return a * a;
}

ConformalMapping::MainForm::MainForm()
{
	InitData();
	InitComponents();
}

ConformalMapping::MainForm::~MainForm()
{
	for (int i = 0; i < m + 1; i++)
	{
		delete[] x[i];
		delete[] y[i];
	}

	delete[] x;
	delete[] y;

	if (components)
	{
		delete components;
	}
}

void ConformalMapping::MainForm::InitData()
{
	x = new double *[m + 1];
	for (int i = 0; i < m + 1; i++)
		x[i] = new double[n + 1];

	y = new double *[m + 1];
	for (int i = 0; i < m + 1; i++)
		y[i] = new double[n + 1];

	Gamma = 0;
	Q = 0;
	k = 0;
}

void ConformalMapping::MainForm::InitialMapping()
{
	double gammaArr[m][n] = {};

	// Initial approximations of boundary nodes
	// Optimize conditions of loops
	for (int j = 0; j < n + 1; j++)
	{
		x[0][j] = (L / (n + 1)) * j;						//	AB
		y[0][j] = 0;										//

		y[m][j] = (-R / (n + 1)) * j - H;					//
		x[m][j] = sqrt(R*R - pow((y[m][j] + H + R), 2));	//	DC	
	}
	for (int i = 0; i < m + 1; i++)
	{
		x[i][0] = 0;										//	AD
		y[i][0] = (-H / (m + 1)) * i;						//

		x[i][n] = ((L - R) / (m + 1)) * i + R;				//	B*C
		y[i][n] = -H - R;									//

		x[i][n] = L;										//	BB*
		y[i][n] = ((-H - R) / (m + 1)) * i;					//
	}

	// Initial approximations of internal nodes
	for (int i = 1; i < m; i++)
	{
		for (int j = 1; j < n; j++)
		{
			x[i][j] = (x[0][j] + x[m][j] + x[i][0] + x[i][n]) / 4.0f;
			y[i][j] = (y[0][j] + y[m][j] + y[i][0] + y[i][n]) / 4.0f;
		}
	}

	double summGammaArr = 0;
	//Initial approximations of GAMMA
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			gammaArr[i][j] = (sqrt(square(x[i + 1][j] - x[i][j]) + square(y[i + 1][j] - y[i][j])) +
				sqrt(square(x[i + 1][j + 1] - x[i][j + 1]) + square(y[i + 1][j + 1] - y[i][j + 1]))) /
				(sqrt(square(x[i][j + 1] - x[i][j]) + square(y[i][j + 1] - y[i][j])) +
					sqrt(square(x[i + 1][j + 1] - x[i + 1][j]) + square(y[i + 1][j + 1] - y[i + 1][j])));

			summGammaArr += gammaArr[i][j];
		}
	}

	Gamma = (1.0f * summGammaArr) / ((m + 1)*(n + 1));
	Q = deltaPhi * (n + 1) / Gamma;

	return;
}

void ConformalMapping::MainForm::Mapping()
{
	k++;

	double gammaArr[m][n] = {};

	// New approximations with Zeydel method

	std::vector<std::vector<double>> tx;
	tx.resize(m + 1, std::vector<double>(n + 1));

	std::vector<std::vector<double>> ty;
	ty.resize(m + 1, std::vector<double>(n + 1));

	for (int i = 0; i < m + 1; i++)
		for (int j = 0; j < n + 1; j++)
		{
			tx[i][j] = x[i][j];
			ty[i][j] = y[i][j];
		}

	for (int i = 1; i < m; i++)
	{
		for (int j = 1; j < n; j++)
		{
			tx[i][j] = (sigma*(x[i + 1][j + 1] - 2 * x[i][j + 1] + x[i - 1][j + 1]) + (1 - 2 * sigma)*(x[i + 1][j] + tx[i - 1][j]) +
				sigma * (tx[i + 1][j - 1] - 2 * tx[i][j - 1] + tx[i - 1][j - 1]) + square(Gamma)*(sigma*(x[i + 1][j + 1] - 2 * x[i + 1][j] + tx[i + 1][j - 1]) +
				(1 - 2 * sigma)*(x[i][j + 1] + tx[i][j - 1]) + sigma * (x[i - 1][j + 1] - 2 * tx[i - 1][j] + tx[i - 1][j - 1]))) /
					(2 * (1 - 2 * sigma)*(1 + square(Gamma)));

			ty[i][j] = (sigma*(y[i + 1][j + 1] - 2 * y[i][j + 1] + y[i - 1][j + 1]) + (1 - 2 * sigma)*(y[i + 1][j] + ty[i - 1][j]) +
				sigma * (ty[i + 1][j - 1] - 2 * ty[i][j - 1] + ty[i - 1][j - 1]) + square(Gamma)*(sigma*(y[i + 1][j + 1] - 2 * y[i + 1][j] + ty[i + 1][j - 1]) +
				(1 - 2 * sigma)*(y[i][j + 1] + ty[i][j - 1]) + sigma * (y[i - 1][j + 1] - 2 * ty[i - 1][j] + ty[i - 1][j - 1]))) /
					(2 * (1 - 2 * sigma)*(1 + square(Gamma)));
		}
	}

	copyBuffers(tx, ty);

	double summGammaArr = 0;
	double summDArr = 0;
	// Finding new values of Gamma, Q, D
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			gammaArr[i][j] = (sqrt(square(x[i + 1][j] - x[i][j]) + square(y[i + 1][j] - y[i][j])) +
				sqrt(square(x[i + 1][j + 1] - x[i][j + 1]) + square(y[i + 1][j + 1] - y[i][j + 1]))) /
				(sqrt(square(x[i][j + 1] - x[i][j]) + square(y[i][j + 1] - y[i][j])) +
					sqrt(square(x[i + 1][j + 1] - x[i + 1][j]) + square(y[i + 1][j + 1] - y[i + 1][j])));

			summGammaArr += gammaArr[i][j];

			summDArr += (sqrt(square(x[i + 1][j + 1] - x[i][j]) + square(y[i + 1][j + 1] - y[i][j]))) /
				(sqrt(square(x[i][j + 1] - x[i + 1][j]) + square(y[i][j + 1] - y[i + 1][j])));
		}
	}

	Gamma = (1.0f * summGammaArr) / ((m + 1)*(n + 1));
	Q = deltaPhi * (n + 1) / Gamma;
	D = summDArr / ((m + 1)*(n + 1));
}

System::Void ConformalMapping::MainForm::drawConformalMapping(System::Drawing::Color c)
{
	array<Point> ^points;
	points = gcnew array<Point>((n + 1)*(m + 1));
	array<Rectangle>^ rect;
	rect = gcnew array<Rectangle>((n + 1)*(m + 1));

	for (int i = 0; i < m + 1; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			points[(n + 1) * i + j] = Drawing::Point((x[i][j] * 22) + 100, (y[i][j] * 22) + 350);
			rect[(n + 1) * i + j] = Drawing::Rectangle(points[(n + 1) * i + j], Drawing::Size(1.0f, 1.0f));
		}
	}
	Graphics ^g = pictureBox1->CreateGraphics();
	g->Clear(BackColor);
	g->DrawRectangles(gcnew Pen(c, 2.0f), rect);
}

System::Void ConformalMapping::MainForm::showAddInfo()
{
	textBox1->Text = Gamma.ToString();
	textBox2->Text = Q.ToString();
	textBox3->Text = k.ToString();
	textBox4->Text = D.ToString();
}

void ConformalMapping::MainForm::InitComponents()
{
	this->pictureBox1 = (gcnew System::Windows::Forms::PictureBox());
	this->button1 = (gcnew System::Windows::Forms::Button());
	this->button2 = (gcnew System::Windows::Forms::Button());
	this->label1 = (gcnew System::Windows::Forms::Label());
	this->label2 = (gcnew System::Windows::Forms::Label());
	this->label3 = (gcnew System::Windows::Forms::Label());
	this->textBox1 = (gcnew System::Windows::Forms::TextBox());
	this->textBox2 = (gcnew System::Windows::Forms::TextBox());
	this->textBox3 = (gcnew System::Windows::Forms::TextBox());
	this->label4 = (gcnew System::Windows::Forms::Label());
	this->textBox4 = (gcnew System::Windows::Forms::TextBox());
	this->button3 = (gcnew System::Windows::Forms::Button());
	(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1))->BeginInit();
	this->SuspendLayout();
	// 
	// pictureBox1
	// 
	this->pictureBox1->BorderStyle = System::Windows::Forms::BorderStyle::FixedSingle;
	this->pictureBox1->Location = System::Drawing::Point(12, 42);
	this->pictureBox1->Name = L"pictureBox1";
	this->pictureBox1->Size = System::Drawing::Size(476, 376);
	this->pictureBox1->TabIndex = 0;
	this->pictureBox1->TabStop = false;
	// 
	// button1
	// 
	this->button1->Location = System::Drawing::Point(12, 12);
	this->button1->Name = L"button1";
	this->button1->Size = System::Drawing::Size(164, 23);
	this->button1->TabIndex = 1;
	this->button1->Text = L"First Aproximation";
	this->button1->UseVisualStyleBackColor = true;
	this->button1->Click += gcnew System::EventHandler(this, &ConformalMapping::MainForm::button1_Click);
	// 
	// button2
	// 
	this->button2->Enabled = false;
	this->button2->Location = System::Drawing::Point(184, 12);
	this->button2->Name = L"button2";
	this->button2->Size = System::Drawing::Size(149, 23);
	this->button2->TabIndex = 2;
	this->button2->Text = L"Next iteration";
	this->button2->UseVisualStyleBackColor = true;
	this->button2->Click += gcnew System::EventHandler(this, &ConformalMapping::MainForm::button2_Click);
	// 
	// label1
	// 
	this->label1->AutoSize = true;
	this->label1->Font = (gcnew System::Drawing::Font(L"Lucida Bright", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
		static_cast<System::Byte>(0)));
	this->label1->Location = System::Drawing::Point(20, 429);
	this->label1->Name = L"label1";
	this->label1->Size = System::Drawing::Size(68, 20);
	this->label1->TabIndex = 3;
	this->label1->Text = L"Gamma:";
	// 
	// label2
	// 
	this->label2->AutoSize = true;
	this->label2->Font = (gcnew System::Drawing::Font(L"Lucida Bright", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
		static_cast<System::Byte>(0)));
	this->label2->Location = System::Drawing::Point(178, 429);
	this->label2->Name = L"label2";
	this->label2->Size = System::Drawing::Size(25, 20);
	this->label2->TabIndex = 4;
	this->label2->Text = L"Q:";
	// 
	// label3
	// 
	this->label3->AutoSize = true;
	this->label3->Font = (gcnew System::Drawing::Font(L"Lucida Bright", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
		static_cast<System::Byte>(0)));
	this->label3->Location = System::Drawing::Point(415, 429);
	this->label3->Name = L"label3";
	this->label3->Size = System::Drawing::Size(22, 20);
	this->label3->TabIndex = 5;
	this->label3->Text = L"k:";
	// 
	// textBox1
	// 
	this->textBox1->BackColor = System::Drawing::SystemColors::Menu;
	this->textBox1->Location = System::Drawing::Point(84, 429);
	this->textBox1->Name = L"textBox1";
	this->textBox1->Size = System::Drawing::Size(88, 20);
	this->textBox1->TabIndex = 6;
	// 
	// textBox2
	// 
	this->textBox2->BackColor = System::Drawing::SystemColors::Menu;
	this->textBox2->Location = System::Drawing::Point(200, 429);
	this->textBox2->Name = L"textBox2";
	this->textBox2->Size = System::Drawing::Size(88, 20);
	this->textBox2->TabIndex = 7;
	// 
	// textBox3
	// 
	this->textBox3->BackColor = System::Drawing::SystemColors::Menu;
	this->textBox3->Location = System::Drawing::Point(434, 429);
	this->textBox3->Name = L"textBox3";
	this->textBox3->Size = System::Drawing::Size(40, 20);
	this->textBox3->TabIndex = 8;
	// 
	// label4
	// 
	this->label4->AutoSize = true;
	this->label4->Font = (gcnew System::Drawing::Font(L"Lucida Bright", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
		static_cast<System::Byte>(0)));
	this->label4->Location = System::Drawing::Point(294, 429);
	this->label4->Name = L"label4";
	this->label4->Size = System::Drawing::Size(25, 20);
	this->label4->TabIndex = 9;
	this->label4->Text = L"D:";
	// 
	// textBox4
	// 
	this->textBox4->BackColor = System::Drawing::SystemColors::Menu;
	this->textBox4->Location = System::Drawing::Point(316, 429);
	this->textBox4->Name = L"textBox4";
	this->textBox4->Size = System::Drawing::Size(88, 20);
	this->textBox4->TabIndex = 10;
	// 
	// button3
	// 
	this->button3->Enabled = false;
	this->button3->Location = System::Drawing::Point(339, 12);
	this->button3->Name = L"button3";
	this->button3->Size = System::Drawing::Size(149, 23);
	this->button3->TabIndex = 11;
	this->button3->Text = L"Print nodes values into file";
	this->button3->UseVisualStyleBackColor = true;
	this->button3->Click += gcnew System::EventHandler(this, &ConformalMapping::MainForm::button3_Click);
	this->button3->Enabled = false;
	// 
	// MainForm
	// 
	this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
	this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
	this->ClientSize = System::Drawing::Size(500, 458);
	this->Controls->Add(this->button3);
	this->Controls->Add(this->textBox4);
	this->Controls->Add(this->label4);
	this->Controls->Add(this->textBox3);
	this->Controls->Add(this->textBox2);
	this->Controls->Add(this->textBox1);
	this->Controls->Add(this->label3);
	this->Controls->Add(this->label2);
	this->Controls->Add(this->label1);
	this->Controls->Add(this->button2);
	this->Controls->Add(this->button1);
	this->Controls->Add(this->pictureBox1);
	this->Name = L"MainForm";
	this->Text = L"MainForm";
	(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1))->EndInit();
	this->ResumeLayout(false);
	this->PerformLayout();

}

System::Void ConformalMapping::MainForm::button1_Click(System::Object^  sender, System::EventArgs^  e)
{
	InitialMapping();
	drawConformalMapping(Color::Green);
	showAddInfo();

	button1->Enabled = false;
	button2->Enabled = true;
	button3->Enabled = true;
}

System::Void ConformalMapping::MainForm::button2_Click(System::Object^  sender, System::EventArgs^  e)
{
	Mapping();
	drawConformalMapping(Color::DarkViolet);
	showAddInfo();
}

System::Void ConformalMapping::MainForm::button3_Click(System::Object^  sender, System::EventArgs^  e)
{
	std::ofstream out("c.txt");

	for (int i = 0; i < m + 1; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			out << std::setw(13) << std::fixed << std::setprecision(4) << "(" << x[i][j] << "," << y[i][j] << ")";
		}
		out << "," << std::endl;
	}
	out << std::endl << std::endl;
}

void ConformalMapping::MainForm::copyBuffers(std::vector<std::vector<double>> &tx, std::vector<std::vector<double>> &ty)
{
	for (int i = 0; i < m + 1; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			x[i][j] = tx[i][j];
			y[i][j] = ty[i][j];
		}
	}

	return;
}

