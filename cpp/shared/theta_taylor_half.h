/*
 * Copyright (c) 2012-2014, Gabriel Leventhal, ETH Zurich
 * All rights reserved.
 *
 * Original MATLAB implementation:
 * Copyright (c) 2010, Nick Higham and Awad Al-Mohy
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer
 *     in the documentation and/or other materials provided with the
 *     distribution.
 *   * Neither the name of the ETH Zurich nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __THETA_TAYLOR__
#define __THETA_TAYLOR__

static double theta[100] = {
  1.95058476281490390628892139801564553636126220e-03,
  7.44366488703259721049221298017073422670364380e-02,
  2.66455438985590931100233547113020904362201691e-01,
  5.24204856864841195474014057253953069448471069e-01,
  8.10268514288732699668571513029746711254119873e-01,
  1.10822540384214773823146060749422758817672729e+00,
  1.41081512609780590494779062282759696245193481e+00,
  1.71493353771263135065794358524726703763008118e+00,
  2.01902706794889441610507674340624362230300903e+00,
  2.32246886135569985398774406348820775747299194e+00,
  2.62492134557707146313987323082983493804931641e+00,
  2.92629977789356443551582742657046765089035034e+00,
  3.22657300191010243395339784910902380943298340e+00,
  3.52578370683355757364552118815481662750244141e+00,
  3.82397968085245620883938499900978058576583862e+00,
  4.12122835471819115582547965459525585174560547e+00,
  4.41759237158269613132688391488045454025268555e+00,
  4.71313685841031748680052260169759392738342285e+00,
  5.00792084833474326899249717826023697853088379e+00,
  5.30200077643247436043338893796317279338836670e+00,
  5.59542758397874262499271935666911303997039795e+00,
  5.88824850385816667142080405028536915779113770e+00,
  6.18050621001630595685583102749660611152648926e+00,
  6.47223984123264983736589783802628517150878906e+00,
  6.76348485176363745097205537604168057441711426e+00,
  7.05427364278112012385690832161344587802886963e+00,
  7.34463569778766167672756637330166995525360107e+00,
  7.63459798213797036225969350198283791542053223e+00,
  7.92418509850397700944313328363932669162750244e+00,
  8.21341959917526232004547637188807129859924316e+00,
  8.50232230999338156607336713932454586029052734e+00,
  8.79091189654228166716620762599632143974304199e+00,
  9.07920599508559611479086015606299042701721191e+00,
  9.36722072360105073585145873948931694030761719e+00,
  9.65497097994831143807914486387744545936584473e+00,
  9.94247064386005519054378964938223361968994141e+00,
  1.02297327700613198686596660991199314594268799e+01,
  1.05167695629449351457651573582552373409271240e+01,
  1.08035919244118723270275950198993086814880371e+01,
  1.10902093188890518149491981603205204010009766e+01,
  1.13766305965402860778112881234847009181976318e+01,
  1.16628646913169795595877076266333460807800293e+01,
  1.19489203001337838117024148232303559780120850e+01,
  1.22348054433121724571265076519921422004699707e+01,
  1.25205274824999897020916250767186284065246582e+01,
  1.28060932750406397673259561997838318347930908e+01,
  1.30915092215469321956788917304947972297668457e+01,
  1.33767814472042481099833821645006537437438965e+01,
  1.36619156356625381931735319085419178009033203e+01,
  1.39469172721236738254901865730062127113342285e+01,
  1.42317919280488371924775492516346275806427002e+01,
  1.45165462025228624298733848263509571552276611e+01,
  1.48011898763128009193223988404497504234313965e+01,
  1.50857392837010362285354858613573014736175537e+01,
  1.53702196120644476451388982241041958332061768e+01,
  1.56546608892378227295694159693084657192230225e+01,
  1.59390837779599063850355378235690295696258545e+01,
  1.62234810794558867996784101705998182296752930e+01,
  1.65078117631577470092452131211757659912109375e+01,
  1.67920195108520111659800022607669234275817871e+01,
  1.70760646886120319720703264465555548667907715e+01,
  1.73599429227210286796889704419299960136413574e+01,
  1.76436773003271305526595824630931019783020020e+01,
  1.79272974337309740633372712181881070137023926e+01,
  1.82108250605575463509921974036842584609985352e+01,
  1.84942718659514326873249956406652927398681641e+01,
  1.87776433519142642580845858901739120483398438e+01,
  1.90609425039085316200271336128935217857360840e+01,
  1.93441715448101732022223586682230234146118164e+01,
  1.96273324788417333763845817884430289268493652e+01,
  1.99104272097071479663554782746359705924987793e+01,
  2.01934575614844895596888818545266985893249512e+01,
  2.04764252844091885208399617113173007965087891e+01,
  2.07593320588701040207979531260207295417785645e+01,
  2.10421794990532440294828120386227965354919434e+01,
  2.13249691563489811585441202623769640922546387e+01,
  2.16077025225443293265925603918731212615966797e+01,
  2.18903810328162968801279930630698800086975098e+01,
  2.21730060685409000598156126216053962707519531e+01,
  2.24555789599312163318245438858866691589355469e+01,
  2.27381009885167451045617781346663832664489746e+01,
  2.30205733894752846424580638995394110679626465e+01,
  2.33029973538277204170299228280782699584960938e+01,
  2.35853740305052106407401879550889134407043457e+01,
  2.38677045385970920676754758460447192192077637e+01,
  2.41499899440742353817768162116408348083496094e+01,
  2.44322312747082683870303299045190215110778809e+01,
  2.47144295295183518135218037059530615806579590e+01,
  2.49965856741739571589278057217597961425781250e+01,
  2.52787006425690847777332237455993890762329102e+01,
  2.55607753385361000653119845082983374595642090e+01,
  2.58428106381739652874784951563924551010131836e+01,
  2.61248073940512384183421090710908174514770508e+01,
  2.64067664445618639490476198261603713035583496e+01,
  2.66886886364635884660856390837579965591430664e+01,
  2.69705748791093711247413011733442544937133789e+01,
  2.72524262703959401221709413221105933189392090e+01,
  2.75342443752984671334615995874628424644470215e+01,
  2.78160318088722284812774887541308999061584473e+01,
  2.80977933868767628666773816803470253944396973e+01,
};

#endif